import flowws
from flowws import Argument as Arg
import freud
import numpy as np
import scipy as sp, scipy.sparse
import sklearn, sklearn.cluster
import plato
import plato.draw as draw

def regularize_box(positions, box_matrix, dtype=None, dimensions=3):
    """ Convert box into a right-handed coordinate frame with
    only upper triangular entries. Also convert corresponding
    positions. Taken from garnett."""
    # First use QR decomposition to compute the new basis
    Q, R = np.linalg.qr(box_matrix)
    Q = Q.astype(dtype)
    R = R.astype(dtype)

    if not np.allclose(Q[:dimensions, :dimensions], np.eye(dimensions)):
        # If Q is not the identity matrix, then we will be
        # changing data, so we have to copy. This only causes
        # actual failures for non-writeable GSD frames, but could
        # cause unexpected data corruption for other cases
        positions = np.copy(positions)

        # Since we'll be performing a quaternion operation,
        # we have to ensure that Q is a pure rotation
        sign = np.linalg.det(Q)
        Q = Q*sign
        R = R*sign

        # First rotate positions and velocities; since they are
        # vectors, we can use the matrix directly. Conveniently,
        # instead of transposing Q we can just reverse the order
        # of multiplication here
        positions = positions.dot(Q)

        # Now we have to ensure that the box is right-handed. We
        # do this as a second step to avoid introducing reflections
        # into the rotation matrix before making the quaternion
        signs = np.diag(np.diag(np.where(R < 0, -np.ones(R.shape), np.ones(R.shape))))
        box = R.dot(signs)
        positions = positions.dot(signs)
    else:
        box = box_matrix

    # Construct the box
    Lx, Ly, Lz = np.diag(box).flatten().tolist()
    xy = box[0, 1]/Ly
    xz = box[0, 2]/Lz
    yz = box[1, 2]/Lz
    box = (Lx, Ly, Lz, xy, xz, yz)
    return positions, box

@flowws.add_stage_arguments
class Projection(flowws.Stage):
    """Project the system using a previously-calculated set of basis vectors.

    This stage wraps the system into the box specified by a given box
    matrix and clusters the particles using sklearn. Either a density
    map or the final cluster centers can be visualized.
    """

    ARGS = [
        Arg('color_scale', None, float, 8,
            help='Scale to adjust colors by in the reduced unit cell density map'),
        Arg('minimum_distance', None, float, .2,
            help='Distance scale to use to cluster projected unit cell particles'),
        Arg('density', '-d', bool, False,
            help='If True, plot a density map instead of clustered points'),
    ]

    def run(self, scope, storage):
        """Project and cluster the particles."""
        box_mat = scope['basis_vectors']

        # if we haven't yet selected a full set of basis vectors, just
        # skip this stage
        if np.abs(np.linalg.det(box_mat)) < 1e-4:
            return

        (positions, box) = regularize_box(scope['position'], box_mat)
        freud_box = freud.box.Box.from_box(box)
        coordinates = plato.math.make_fractions(box, positions)
        coordinates %= 1.
        positions = self.projected_positions = \
            plato.math.fractions_to_coordinates(box, coordinates)
        types = self.projected_types = scope['type']
        num_types = len(np.unique(types))

        nq = freud.locality.AABBQuery(freud_box, positions)
        qargs = dict(r_max=self.arguments['minimum_distance'], exclude_ii=True)
        nlist = nq.query(positions, qargs).toNeighborList()

        rijs = (positions[nlist.point_indices] -
                positions[nlist.query_point_indices])
        rijs = freud_box.wrap(rijs)
        rs = np.sqrt(np.sum(rijs**2, axis=-1))

        distance_matrix = sp.sparse.coo_matrix(
            (rs, (nlist.query_point_indices, nlist.point_indices)),
            (positions.shape[0], positions.shape[0])).tocsr()
        dbscan = sklearn.cluster.DBSCAN(
            self.arguments['minimum_distance'], metric='precomputed')
        clusters = dbscan.fit_predict(distance_matrix)

        possible_clusters = list(sorted(np.unique(clusters)))
        cluster_indices = [np.where(clusters == idx)[0] for idx in possible_clusters]

        if possible_clusters[0] == -1:
            unclustered_indices = cluster_indices[0]
            cluster_indices = cluster_indices[1:]
        else:
            unclustered_indices = np.array([], dtype=np.uint32)

        cluster_centers = np.zeros((len(cluster_indices), 3), dtype=np.float32)
        cluster_type_fractions = np.zeros((len(cluster_indices), num_types), dtype=np.float32)
        cluster_histogram = np.array([len(ix) for ix in cluster_indices], dtype=np.uint32)

        for (i, indices) in enumerate(cluster_indices):
            ref = positions[indices[0]]
            shifted_positions = positions[indices] - ref[np.newaxis, :]
            shifted_positions = freud_box.wrap(shifted_positions)

            cluster_centers[i] = np.mean(shifted_positions, axis=0) + ref
            cluster_type_fractions[i] = np.bincount(types[indices], minlength=num_types)/len(indices)

        scope['box'] = self.box = box
        scope['position'] = self.positions = cluster_centers
        scope['type'] = self.types = np.argmax(cluster_type_fractions, axis=-1)
        scope.setdefault('visuals', []).append(self)

    def draw_plato(self):
        min_length = min(*self.box[:3])
        pixel_scale = 600/min_length
        size = (min_length*4./3, min_length)

        box_args = dict(zip(['Lx', 'Ly', 'Lz', 'xy', 'xz', 'yz'], self.box))
        box_prim = draw.Box(**box_args)
        box_prim.width = min_length*.05

        if self.arguments['density']:
            box_prim.color = (1, 1, 1, 1)
            color_scale = (self.arguments['color_scale']*len(self.positions)/
                           len(self.projected_positions))

            prim = draw.Spheres(positions=self.projected_positions, diameters=.05)
            colors = np.zeros((len(prim.positions), 4))
            colors[:, :3] = plato.cmap.cubeellipse(self.projected_types.astype(np.float32))
            prim.colors = colors*color_scale

            features = dict(ambient_light=0, directional_light=(0, 0, -1),
                            additive_rendering=True)
            scene = draw.Scene([box_prim, prim], features=features)
        else:
            colors = np.ones((len(self.positions), 4))
            colors[:, :3] = plato.cmap.cubeellipse(self.types.astype(np.float32))

            prim = draw.Spheres(positions=self.positions, colors=colors)
            scene = draw.Scene([box_prim, prim])

        scene.size = size
        scene.pixel_scale = pixel_scale

        return scene
