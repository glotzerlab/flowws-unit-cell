import flowws
from flowws import Argument as Arg
import freud
import numpy as np
import plato
import plato.draw.vispy as draw
import rowan

def circle_patterns(locations, radii, Npoints=128, z=0):
    locations = np.array(locations)

    thetas = np.linspace(0, 2*np.pi, Npoints, endpoint=False)
    circle_template = np.zeros((Npoints, 3))
    circle_template[:, 0] = np.cos(thetas)
    circle_template[:, 1] = np.sin(thetas)

    result = []
    for (location, radius) in zip(locations, radii):
        result.append(location[np.newaxis, :] + circle_template*radius)
        result[-1][:, 2] = z

    return np.concatenate(result, axis=0)

def find_best_orientation(bonds, quat):
    bonds = bonds.copy()
    # rotate bonds onto [0, 0, 1]
    bonds[bonds[:, 2] < 0] *= -1
    mean = np.mean(bonds, axis=0)
    norm_mean = mean/np.sqrt(np.sum(mean**2))

    target = [0, 0, 1]

    axis = np.cross(norm_mean, target)
    axis /= np.sqrt(np.sum(axis**2))
    halfangle = np.arccos(norm_mean[2])/2
    adjustment_quat = np.array([np.cos(halfangle)] + (np.sin(halfangle)*axis).tolist())

    result = plato.math.quatquat(adjustment_quat, quat)
    return result, mean

class FixedSpherePoints(draw.SpherePoints):
    @property
    def rotation(self):
        return None

    @rotation.setter
    def rotation(self, value):
        pass

# make this available for scene.convert()
draw.FixedSpherePoints = FixedSpherePoints

class DirectionVisual:
    def update(self, angle_tolerance, circle_positions, bonds, rotation,
                 intensity=1e4):
        self.angle_tolerance = angle_tolerance
        self.circle_positions = circle_positions
        self.bonds = bonds
        self.rotation = rotation

        self.make_reference_lines()
        self.bond_prim = draw.SpherePoints(
            points=bonds, on_surface=True, intensity=intensity)
        self.circle_prim = draw.SpherePoints(on_surface=False, intensity=2e1)

    def make_reference_lines(self):
        self.R = np.sin(self.angle_tolerance*np.pi/180)
        reference_lines = []

        reference_Rs = np.linspace(self.R, 1, 128)
        for theta in [0, np.pi/4, np.pi/6, np.pi/2]:
            xs = reference_Rs*np.cos(theta)
            ys = reference_Rs*np.sin(theta)

            reference_lines.append(np.array([xs, ys, np.zeros_like(xs)], dtype=np.float32).T)
            reference_lines.append(np.array([-xs, ys, np.zeros_like(xs)], dtype=np.float32).T)
            reference_lines.append(np.array([xs, -ys, np.zeros_like(xs)], dtype=np.float32).T)
            reference_lines.append(np.array([-xs, -ys, np.zeros_like(xs)], dtype=np.float32).T)
        reference_lines = np.concatenate(reference_lines, axis=0)
        reference_lines = np.ascontiguousarray(reference_lines, dtype=np.float32)

        self.reference_prim = FixedSpherePoints(
            points=reference_lines, on_surface=False, intensity=1e2)

    def draw_plato(self):
        direction_scene = draw.Scene(
            [self.bond_prim, self.circle_prim, self.reference_prim], zoom=14,
            features=dict(additive_rendering=dict(invert=True)),
            rotation=self.rotation)

        circle_positions = (np.concatenate(self.circle_positions) if self.circle_positions
                            else np.array([0, 0, 0]))
        self.circle_prim.points = circle_positions

        return direction_scene

@flowws.add_stage_arguments
class BasisSelection(flowws.Stage):
    """Select directions and distances to form the basis vectors for a unit cell.

    This stage produces two visuals: a bond orientational order
    diagram that can be used to select symmetric directions, and a
    "cylindrical RDF" that measures bonds along each given
    direction. Together, these can be used to select the direction and
    length of the three basis vectors for the unit cell.
    """
    ARGS = [
        Arg('orientations', '-d', [(float, float, float, float)], [],
            help='Quaternions specifying orientations for basis vectors (in the (0, 0, 1) direction)'),
        Arg('r_max', '-r', float, 3,
            help='Maximum distance to consider for binds to select for the unit cell basis vectors'),
        Arg('angle_tolerance', '-a', float, 5,
            help='Angle tolerance for selecting bonds (in degrees)'),
        Arg('rdf_bins', None, int, 128,
            help='Number of bins to use for cylindrical RDF'),
        Arg('x_direction', '-x', int, 0,
            help='Candidate direction to take as the x direction in the final unit cell'),
        Arg('x_min', None, float, 0,
            help='Minimum distance to take bonds from the cylindrical RDF for the x basis vector'),
        Arg('x_max', None, float, 1,
            help='Maximum distance to take bonds from the cylindrical RDF for the x basis vector'),
        Arg('y_direction', '-y', int, 1,
            help='Candidate direction to take as the y direction in the final unit cell'),
        Arg('y_min', None, float, 0,
            help='Minimum distance to take bonds from the cylindrical RDF for the y basis vector'),
        Arg('y_max', None, float, 1,
            help='Maximum distance to take bonds from the cylindrical RDF for the y basis vector'),
        Arg('z_direction', '-z', int, 2,
            help='Candidate direction to take as the z direction in the final unit cell'),
        Arg('z_min', None, float, 0,
            help='Minimum distance to take bonds from the cylindrical RDF for the z basis vector'),
        Arg('z_max', None, float, 1,
            help='Maximum distance to take bonds from the cylindrical RDF for the z basis vector'),
    ]

    def __init__(self, *args, **kwargs):
        self.direction_visual = DirectionVisual()
        super().__init__(*args, **kwargs)

    def run(self, scope, storage):
        """Display the interactive direction- and distance-selection visuals."""
        positions = scope['position']
        types = scope['type']

        box = freud.box.Box.from_box(scope['box'])
        aq = freud.AABBQuery(box, positions)
        args = dict(exclude_ii=True, mode='ball', r_max=self.arguments['r_max'])
        nlist = aq.query(positions, args).toNeighborList()
        bonds = positions[nlist.point_indices] - positions[nlist.query_point_indices]
        bonds = box.wrap(bonds)

        angle_tol = self.arguments['angle_tolerance']*np.pi/180
        R = np.sin(angle_tol)

        normalized_bonds = bonds/nlist.distances[:, np.newaxis]
        orientations = np.array(self.arguments['orientations'])
        orientations /= np.linalg.norm(orientations, axis=-1, keepdims=True)

        distances = np.linalg.norm(bonds, axis=-1)

        basis = np.eye(3)

        bond_filter = np.ones((bonds.shape[0],), dtype=np.bool)
        target_quat = np.array((1., 0, 0, 0))

        optimized_orientations = []
        found_vectors = []
        found_bonds = []
        circle_positions = []
        for q in orientations:
            qconj = plato.math.quatconj(q)

            rotated_bonds = plato.math.quatrot(q[np.newaxis, :], bonds)
            filt_thetas = np.arctan2(np.sqrt(rotated_bonds[:, 0]**2 + rotated_bonds[:, 1]**2),
                                     rotated_bonds[:, 2])
            filt_thetas[filt_thetas > np.pi/2] -= np.pi
            filt = np.abs(filt_thetas) < angle_tol

            filtered_bonds = rotated_bonds[filt]

            found_bonds.append(filtered_bonds)

            if filtered_bonds.shape[0]:
                (target_quat, mean) = find_best_orientation(filtered_bonds, q)
            else:
                target_quat = q

            qconj = plato.math.quatconj(target_quat)

            found_vectors.append(plato.math.quatrot(qconj, [0, 0, 1]))
            optimized_orientations.append(target_quat)

            circle = circle_patterns([(0, 0, 0)], [R], z=-1)
            circle = plato.math.quatrot(qconj, circle)
            circle_positions.append(circle)

            bond_filter *= np.logical_not(filt)

        self.direction_radii = [np.linalg.norm(b, axis=-1) for b in found_bonds]

        for i, key in enumerate('xyz'):
            direction_name = key + '_direction'
            min_name, max_name = key + '_min', key + '_max'

            direction_index = self.arguments[direction_name]
            r_min, r_max = self.arguments[min_name], self.arguments[max_name]
            r_min, r_max = sorted([r_min, r_max])

            try:
                direction = found_vectors[direction_index]
            except IndexError:
                break

            candidate_distances = self.direction_radii[direction_index]
            filt = np.logical_and(
                candidate_distances >= r_min, candidate_distances < r_max)
            if np.any(filt):
                candidate_distances = candidate_distances[filt]
            distance = np.mean(candidate_distances)

            basis[:, i] = direction*distance

        if np.linalg.det(basis) < 0:
            basis[:, 0] *= -1

        self.direction_visual.update(
            self.arguments['angle_tolerance'], circle_positions, bonds[bond_filter],
            target_quat.copy())

        self.basis = scope['basis_vectors'] = basis
        scope.setdefault('visuals', []).append(self.direction_visual)
        scope['visuals'].append(self)

        self.gui_actions = [
            ('Select bonds', self._select_current_position),
            ('Undo last selection', self._undo_selection),
        ]

    def draw_matplotlib(self, figure):
        import matplotlib

        ax = figure.add_subplot()
        colors = matplotlib.rcParams['axes.prop_cycle'].by_key()['color']

        for i, (rs, color) in enumerate(zip(self.direction_radii, colors)):
            (hist, _) = np.histogram(rs, bins=self.arguments['rdf_bins'])

            for key in 'xyz':
                direction_name = key + '_direction'
                if self.arguments[direction_name] != i:
                    continue
                min_name, max_name = key + '_min', key + '_max'
                r_min, r_max = self.arguments[min_name], self.arguments[max_name]
                r_min, r_max = sorted([r_min, r_max])

                r_taken = np.linalg.norm(self.basis[:, 'xyz'.index(key)])
                ax.vlines([r_min, r_max], 0, np.max(hist)*1.05,
                          linestyles='dashed', color=color)
                ax.vlines([r_taken], 0, np.max(hist)*1.05, linestyles='solid',
                          color=color)

            ax.hist(
                rs, bins=self.arguments['rdf_bins'], alpha=.5,
                label='Direction {}'.format(i), color=color)

        ax.legend()
        ax.set_xlabel('$R$')
        ax.set_ylabel('$Count(R)$')

    def _select_current_position(self, scope, storage):
        plato_scene = scope['visual_objects'][self.direction_visual]

        self.arguments['orientations'].append(plato_scene.rotation.copy())

        if scope.get('rerun_callback', None) is not None:
            scope['rerun_callback']()

    def _undo_selection(self, scope, storage):
        if self.arguments['orientations']:
            self.arguments['orientations'].pop()

        if scope.get('rerun_callback', None) is not None:
            scope['rerun_callback']()
