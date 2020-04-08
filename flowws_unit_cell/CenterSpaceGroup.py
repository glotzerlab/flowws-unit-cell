from collections import Counter
import flowws
from flowws import Argument as Arg
import numpy as np
import plato
import spglib

def identify_spacegroup(spg_cell, max_iterations=200, minimum_distance=.9):
    """This function aims to identify the best spacegroup
    based on allowing some amount of deviation. Accepts a spglib
    compatible cell tuple as the argument."""

    # Note that while precision to spglib is in cartesian distance,
    # input positions are fractional coordinates
    prec = 5*minimum_distance
    precs = []
    grps = []
    grp = 'None'
    highest_symmetry_group = 'None'
    max_group = 0

    # Try a range of precisions and record the determined spacegroups
    counter = 0
    while grp is None or grp.split()[-1] not in ['(1)', '(2)']:
        counter += 1
        if counter > max_iterations:
            break
        grp = spglib.get_spacegroup(spg_cell, symprec=prec)
        grps.append(grp)
        precs.append(prec)
        prec /= 2

        if grp is not None:
            group_num = int(grp.split()[-1].replace('(', '').replace(')', ''))
            if group_num > max_group:
                max_group = group_num
                highest_symmetry_group = grp

    if all(g is None for g in grps):
        raise ValueError("No symmetry groups found!")
    # One measure of best group is the highest symmetry group
    highest_symmetry_group_prec = precs[::-1][grps[::-1].index(highest_symmetry_group)]

    # An alternative measure is the most commonly occurring group
    counts = Counter(grps)
    if None in counts:
        del counts[None]
    most_common_group = counts.most_common(1)[0][0]
    most_common_group_prec = precs[::-1][grps[::-1].index(most_common_group)]

    return {'common': (most_common_group, most_common_group_prec),
            'highest': (highest_symmetry_group, highest_symmetry_group_prec),
            'histogram': counts}

def matrix_to_box(bm):
    """Convert a box matrix into a box object"""
    bm = np.asarray(bm)
    Lx, Ly, Lz = np.diag(bm).flatten().tolist()
    xy = bm[0, 1]/Ly
    xz = bm[0, 2]/Lz
    yz = bm[1, 2]/Lz
    box = (Lx, Ly, Lz, xy, xz, yz)

    return box

def standardize_cell(spg_cell, best_prec):
    """Convert cell into its standard crystallographic representation"""
    dataset = spglib.get_symmetry_dataset(spg_cell, best_prec)
    box = matrix_to_box(dataset['std_lattice'].T)
    positions = dataset['std_positions']

    return (box, positions, dataset['std_types'])

@flowws.add_stage_arguments
class CenterSpaceGroup(flowws.Stage):
    """Attempt to automatically detect the space group of the system and center it accordingly."""
    ARGS = [
        Arg('minimum_distance', '-d', float, .9,
            help='Precision with which to merge points transformed by space group transformations'),
        Arg('use_types', '-t', bool, True,
            help='Use type information when centering the unit cell'),
    ]

    def run(self, scope, storage):
        """Detect the space group and center the system."""
        box = scope['box']
        boxmat = plato.math.box_to_matrix(box)
        fractions = plato.math.make_fractions(box, scope['position'])
        types = scope['type']

        if not self.arguments['use_types']:
            types = np.zeros_like(types)

        spglib_cell = (boxmat, fractions, types)
        try:
            self.spg_info = identify_spacegroup(
                spglib_cell, max_iterations=64,
                minimum_distance=self.arguments['minimum_distance'])
        except ValueError:
            # didn't work; don't display things
            return

        (box, fractions, types) = standardize_cell(
            spglib_cell, self.spg_info['common'][1])

        scope['box'] = box
        scope['position'] = plato.math.fractions_to_coordinates(box, fractions)
        scope['type'] = types
        scope.setdefault('visuals', []).append(self)

    def draw_matplotlib(self, figure):
        ax = figure.add_subplot()

        histogram = self.spg_info['histogram']
        keys = list(histogram)
        counts = [histogram[k] for k in keys]
        xs = np.arange(len(keys))
        ax.bar(xs, counts)

        ax.set_title('Space group prevalence')
        ax.set_xticks(xs)
        ax.set_xticklabels(keys)
        ax.set_ylabel('Count')
