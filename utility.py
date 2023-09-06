import open3d as o3d
import numpy as np
import matplotlib.pyplot as plt
import warnings
from matplotlib import cm
from scipy.sparse import csr_matrix, lil_matrix, diags, eye
from scipy.sparse.linalg import eigs, eigsh
from scipy.linalg import eig
from random import sample


def create_test_mesh():

    theta = np.linspace(0, 2*np.pi, 7)[:-1]
    x, y, z = np.cos(theta), np.zeros(6), np.sin(theta)

    vertices = np.vstack((x,y,z)).T
    vertices = np.vstack((np.array([0,1,0]), vertices))

    triangles = np.array([
        [0,1,2], [0,2,3], [0,3,4], [0,4,5], [0,5,6], [0,6,1],
        [0,2,1], [0,3,2], [0,4,3], [0,5,4], [0,6,5], [0,1,6]
    ])

    delta = delta_coordinates_single(0, vertices, triangles)

    o3d.visualization.draw_geometries([
        o3d.geometry.TriangleMesh(
            o3d.utility.Vector3dVector(vertices),
            o3d.utility.Vector3iVector(triangles)
        ),

        o3d.geometry.LineSet(
            o3d.utility.Vector3dVector(np.array([vertices[0,:], vertices[0,:]-delta])),
            o3d.utility.Vector2iVector(np.array([[0,1]]))
        ).paint_uniform_color(np.array([0,1,0]))
    ])

def adjacency_matrix_dense(triangles, num_vertices=None):

    if num_vertices is None:
        num_vertices = triangles.max()+1
    
    #initializing sparse matrix
    adj_matrix = np.zeros((num_vertices, num_vertices), dtype=np.uint16)

    #iterating triangles to populate the array
    for tri in triangles:

        v1, v2, v3 = tri
        adj_matrix[v1, v2] = 1
        adj_matrix[v2, v1] = 1
        adj_matrix[v2, v3] = 1
        adj_matrix[v3, v2] = 1
        adj_matrix[v3, v1] = 1
        adj_matrix[v1, v3] = 1

    return adj_matrix

def adjacency_matrix_sparse(triangles, num_vertices = None):

    if num_vertices is None:
        num_vertices = triangles.max()+1

    #initializing sparse matrix
    adj_matrix = lil_matrix((num_vertices, num_vertices), dtype=np.uint16)

    #iterating triangles to populate the array
    for tri in triangles:

        v1, v2, v3 = tri
        adj_matrix[v1, v2] = 1
        adj_matrix[v2, v1] = 1
        adj_matrix[v2, v3] = 1
        adj_matrix[v3, v2] = 1
        adj_matrix[v3, v1] = 1
        adj_matrix[v1, v3] = 1

    #converting to csr
    adj_matrix = adj_matrix.tocsr()

    return adj_matrix

def degree_matrix(adj, exponent=1):

    num_vertices = adj.shape[0]
    diagonals = np.zeros(num_vertices)

    if exponent==1:
        for i in range(num_vertices):
            diagonals[i] = adj[i,:].toarray().sum()
        return diags(diagonals, format="csr", dtype=np.int32)
    else:
        for i in range(num_vertices):
            diagonals[i] = adj[i,:].toarray().sum().astype(np.float32)**exponent
        return diags(diagonals, format="csr", dtype=np.float32)

def delta_coordinates_single(idx, vertices, triangles, k=1):

    vi = vertices[idx]

    neighbors = k_ring_adjacency(idx, triangles, k)
    delta = vi - vertices[neighbors,:].mean(0)

    return delta

def delta_coordinates(vertices, triangles, use_laplacian=True):

    if use_laplacian:
        L = random_walk_laplacian(triangles)
        delta = L @ vertices
    else:
        delta = np.zeros_like(vertices)
        for i, vi in enumerate(vertices):
            neighbors = k_ring_adjacency(i, triangles, 1)
            delta[i] = vi - vertices[neighbors, :].mean(0)
    
    return delta

def adjacency_matrix(triangles, num_vertices = None):

    if num_vertices is None:
        num_vertices = triangles.max()+1

    #initializing sparse matrix
    adj_matrix = lil_matrix((num_vertices, num_vertices), dtype=np.uint16)

    #iterating triangles to populate the array
    for tri in triangles:

        v1, v2, v3 = tri
        adj_matrix[v1, v2] = 1
        adj_matrix[v2, v1] = 1
        adj_matrix[v2, v3] = 1
        adj_matrix[v3, v2] = 1
        adj_matrix[v3, v1] = 1
        adj_matrix[v1, v3] = 1

    #converting to csr
    adj_matrix = adj_matrix.tocsr()

    return adj_matrix

def k_ring_recursive(idx, triangles, k = 1):

    #terminating condition, k=0
    if not k:
        return np.array([])

    #from a to np.array([a])
    if isinstance(idx, int):
        idx = np.array([idx], dtype=np.uint16)

    new_ids = np.array([], dtype=np.uint16)

    #iterating all vertex ids
    for id in idx:
        #iterating all triangles
        for t in triangles:
            #if id is found in a triangle, add the other two vertices
            if id in t:
                new_ids = np.hstack((new_ids, t[t - id != 0]))

    #discarding duplicates
    new_ids = np.unique(new_ids)

    return np.unique(np.hstack((new_ids, k_ring_recursive(new_ids, triangles, k-1)))).astype(np.uint32)

def k_ring_adjacency(idx, triangles, k=1, num_vertices=None):

    adj_matrix = adjacency_matrix(triangles, num_vertices)

    #kth power
    adj_matrix = adj_matrix ** k

    #neighbors of specified index are the indices of non-zero elements of that row
    neighbors = adj_matrix[idx, :].toarray()

    return neighbors.nonzero()[1]

def sample_colormap(scalars, name="inferno"):

    avail_maps = ["inferno", "magma", "viridis", "cividis"]

    if name not in avail_maps:
        warnings.warn(f"Only {avail_maps} colormaps are supported. Using inferno.")
        name = "inferno"

    colormap = cm.get_cmap(name, 12)
    colors = colormap(scalars)

    return colors[:,:-1]

def graph_laplacian(triangles):

    num_vertices = triangles.max()+1

    A = adjacency_matrix(triangles, num_vertices=num_vertices)
    D = degree_matrix(A, exponent=1)

    L = D - A

    return L

def random_walk_laplacian(triangles, subtract=True):

    num_vertices = triangles.max()+1

    A = adjacency_matrix(triangles, num_vertices=num_vertices)
    Dinv = degree_matrix(A, exponent=-1)

    if subtract:
        L = eye(num_vertices, num_vertices, 0) - Dinv @ A
    else:
        L = Dinv @ A

    return L

def laplacian_smoothing(triangles, vertices, smoothL):
    delta = delta_coordinates(vertices, triangles, use_laplacian = True)
    vertices -= smoothL * delta
    return vertices

def taubin_smoothing(triangles, vertices, smoothL, smoothM):
    delta = delta_coordinates(vertices, triangles, use_laplacian=True)
    vertices -= smoothL * delta
    delta = delta_coordinates(vertices, triangles, use_laplacian=True)
    vertices += smoothM * delta
    return vertices

def compute_len_stats(vertices:o3d.utility.Vector3dVector):
    pcd = o3d.geometry.PointCloud(vertices)
    d = np.asarray(pcd.compute_nearest_neighbor_distance())
    return d.min(), d.mean(), d.max()

def add_shortcut_edges(triangles, vertices):
    #Something something adjacency matrix, Rotations and other weird stuff?
    for triangle in triangles:
        pass
    pass

def geodesic_dijkstra(vertices, mesh):
    
    g_u = float('inf') * np.ones(len(vertices))
    base_vertices_num = round(len(g_u)/6)
    # base_vertices
    base_vertices = np.asarray(mesh.sample_points_uniformly(base_vertices_num).points)
    #get base indexes (N complexity)
    #

    pass


if __name__ == "__main__":

    create_test_mesh()