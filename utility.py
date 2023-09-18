import open3d as o3d
import numpy as np
import matplotlib.pyplot as plt
import warnings
from matplotlib import cm
from scipy.sparse import csr_matrix, lil_matrix, diags, eye
from scipy.sparse.linalg import eigs, eigsh
from scipy.linalg import eig
from random import sample
import heapq



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



def geodesic_dijkstra(vertices, triangles, S_area):
    V_LIST_NON_EMPTY = 1
    vertlen = len(vertices)
    r_threshold = np.sqrt(0.005 * S_area)
    print("r:", r_threshold)
    A_matrix = adjacency_matrix(triangles)

   
    unvisited_vertices = np.zeros(vertlen) #idxs of visited vertices
    mu_values = np.zeros(vertlen)
    g_values = []
    vlist = []
    # vlist = vertices.copy()
    # vlist = heapq.heapify(vlist)
    #Base0 is vertex 0. Initialize first base to 0
    last_index = 0
    # g_u[last_index] = 0
    # visited[index] = 1
    base_points = []
    
    vlist.append(last_index)
    bases_counter = 0
    base_points_length = 0
    # vlist = heapq.heapify(vlist)
    runtime = True
    while(runtime):
        j = last_index
        temporary_point = None
        while (j < len(unvisited_vertices)):
            if (unvisited_vertices[j]!= 1):
                temp = vertices[j]
                last_index = j
                break
            j+= 1
        if(temporary_point != None):
            temp_int = j
            if(base_points_length >= len(base_points)):
                #weird expansions
                pass
            base_points[base_points_length] = temporary_point
            minHeapVLIST = [[float('inf'), vertex_index] for vertex_index in range(vertlen)]
            g_values[base_points_length] = calculateShortestPath(temp_int, minHeapVLIST, vertlen, A_matrix,
                                                                  r_threshold, unvisited_vertices)
            #g_values is 2D matrix where oneD is bases and second is distance of all vertices
            #to that base.
            



    # while sum(visited) < len(visited):
    #     while V_LIST_NON_EMPTY:
    #         u = vlist.pop()
    #         neighbours_u = A_matrix[u, :].toarray()
    #         neighbours_u = neighbours_u.nonzero()[1]
    #         for ua in neighbours_u:
    #             dist_check = g_u[u] + np.linalg.norm(vertices[u] - vertices[ua])
    #             visited[ua] = 1
    #             if (g_u[ua] > dist_check):
    #                 g_u[ua] = dist_check
    #                 if (g_u[ua] < r_threshold):
    #                     if(len(vlist) == 0 or g_u[vlist[-1]] <= g_u[ua]):
    #                         vlist.append(ua)
    #                     elif (g_u[vlist[-1]] > g_u[ua]):
    #                         vlist.insert(-1, ua)
    #         V_LIST_NON_EMPTY = len(vlist)
    #     new_base = np.where(visited == 0)[0][0]
    #     visited[new_base] = 1
    #     bases_counter += 1
    #     vlist.append(new_base)
    # print("Based on ", bases_counter, " bases.")
    # print("vertices: ", len(vertices))
    # m_u = sum(g_u) 
    ##TBC

def calculateShortestPath(base_vertex_index, VLIST, vertlen,
                           A_matrix, vertices, r_threshold,
                             unvisited_vertices):
    #Initialize g(u) to infinity for all indexes
    KEY = 0
    INDEX = 1

    g_bu = float('inf') * np.ones(vertlen)
    VLIST[base_vertex_index][0] = 0
    heapq.heapify(VLIST)
    while(len(VLIST)!= 0):
        smallest = VLIST.pop()
        g_V = smallest[KEY]
        if(g_bu[smallest[INDEX]] > smallest[KEY]):
            g_bu[smallest[INDEX]] = smallest[KEY]
        neighbours_u = A_matrix[smallest[INDEX], :].toarray()
        neighbours_u = neighbours_u.nonzero()[1]
        for neighbour in neighbours_u:
            length_VVa = np.linalg.norm(vertices[smallest[INDEX]] - vertices[neighbour])
            g_Va = g_bu[neighbour]

        if(g_Va > g_V + length_VVa):
            g_bu[neighbour] = g_V + length_VVa
            #Decrease Key:
            if (VLIST[neighbour][KEY] < g_bu[neighbour]):
                VLIST[neighbour][KEY] = g_bu[neighbour]
                heapq.heapify(VLIST)

    base_ver = vertices[base_vertex_index]
    #create array of vertices in the area

    points_in_area_length = 0
    points_in_area = []
    #Check if vertex is within area bounds
    for vertex_index, distance in enumerate(g_bu):
        if(distance <= r_threshold):
            #check if it has already been included in the area
            if (unvisited_vertices[vertex_index] and 
                vertex_index != base_vertex_index):
                points_in_area[points_in_area_length] = vertex_index
                points_in_area_length+= 1

            unvisited_vertices[vertex_index] = 0
        #Here come some expansions, is this just appends?
        # Can you do those in a more sophisticated manner?        





if __name__ == "__main__":

    create_test_mesh()