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
import numpy_tricks as npt
import time



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
    
    INF = 100000
    vertlen = len(vertices)
    r_threshold = np.sqrt(0.005 * S_area)
    print("r:", r_threshold)
    A_matrix = adjacency_matrix(triangles)
    # base_areas = np.zeros(10)
   
    unvisited_vertices = np.ones(vertlen, dtype = np.int8) #idxs of visited vertices
    base_default = 200

    g_values = -np.ones((base_default, vertlen), dtype = np.float128)
    mu_values = np.empty(vertlen, dtype = np.float128)
    
    

    base_points = -np.ones(base_default, dtype = np.int32)
    base_points_length = 0
    base_areas = np.empty(base_default, dtype = np.int32)
    # vlist = []
    # vlist = vertices.copy()
    # vlist = heapq.heapify(vlist)
    #Base0 is vertex 0. Initialize first base to 0
    last_index = 0
    # g_u[last_index] = 0
    # visited[index] = 1
    
    # debug counters and timers
    debug_counters = [0, 0, 0]
    heap_restructuring_time = 0
    base_area_timer = 0
    shortest_path_timer = 0
    # vlist.append(last_index)
    
    # vlist = heapq.heapify(vlist)
    runtime = True
    while(runtime):
        j = last_index
        temporary_point = None
        while (j < vertlen):
            if (unvisited_vertices[j]):
                temporary_point = vertices[j]
                last_index = j
                break
            j+= 1
        
        if(temporary_point is not None):
            temp_int = j

            if(base_points_length >= (len(base_points) - 1)):
                #weird expansions
                print("You're getting cooked.")
                npt.expand_base_empty(base_points, 1)
                npt.expand_base_empty(base_areas, 1)
                npt.expand_rows(g_values, 1)

            #Contrast to the original implementation, we use the index to represent base. 
            base_points[base_points_length] = temp_int #was temporary_point

            minHeapVLIST = [[INF, vertex_index] for vertex_index in range(vertlen)]
            
            shortest_path_timer_start = time.process_time()
            #populate base point row with its g_values
            g_values[base_points_length, : vertlen], heap_restructuring_time, base_area_timer = calculateShortestPath(temp_int, minHeapVLIST, vertlen, A_matrix,
                                                                  vertices, triangles, r_threshold, unvisited_vertices,
                                                                    base_areas, base_points_length, debug_counters,
                                                                      heap_restructuring_time, base_area_timer)
            shortest_path_timer_end = time.process_time()
            shortest_path_timer += shortest_path_timer_end - shortest_path_timer_start
            #g_values is 2D matrix where oneD is bases and second is distance of all vertices
            #to that base.
            base_points_length += 1
            if((base_points_length % 10) == 0):
                print(base_points_length)
        elif(temporary_point is None):
            runtime = False
    

    ## calc mu values and normalize them
    for value_index in range(len(mu_values)):
        mu_values[value_index] = calculate_mu(g_values, base_points_length, value_index, base_areas)
    
    # normalize_mu(mu_values)
    print("Total number of base points:", base_points_length)

    print("Total Heap restructuring time:", heap_restructuring_time)
    print("Total cpu time for Base Area Calculations: ", base_area_timer)
    print("Total cpu time for shortest path timer: ", shortest_path_timer)




def calculateShortestPath(base_vertex_index, VLIST, vertlen,
                           A_matrix, vertices, triangles, r_threshold,
                             unvisited_vertices, base_areas, base_points_length,
                               debug_counters, heap_restructuring_time, base_area_timer):
    #Initialize g(u) to infinity for all indexes
    KEY = 0
    INDEX = 1
    INF = 1000000
    g_bu = INF * np.ones(vertlen)

    #VLIST is 2D, [[key = value, index = vertex_index]]
    VLIST[base_vertex_index][0] = 0
    heapq.heapify(VLIST)
    while(len(VLIST)):
        smallest = heapq.heappop(VLIST)
        g_V = smallest[KEY]

        if(g_bu[smallest[INDEX]] > smallest[KEY]):
            g_bu[smallest[INDEX]] = smallest[KEY]

        neighbours_u = A_matrix[smallest[INDEX], :].toarray()
        neighbours_u = neighbours_u.nonzero()[1]
        for neighbour in neighbours_u:
            length_VVa = np.linalg.norm(vertices[smallest[INDEX]] - vertices[neighbour])
            g_Va = g_bu[neighbour]
            #Main check
            if(g_Va > g_V + length_VVa):
                g_bu[neighbour] = g_V + length_VVa

            #Decrease Key:
            #This method has linear search accross large heap, 
            #But I believe the converter idea is equally stinky.
            heap_time_start = time.process_time()

            decrease_key(VLIST, neighbour, g_bu[neighbour])
            
            heap_time_end = time.process_time()
            heap_restructuring_time += heap_time_end - heap_time_start

            # if (VLIST[neighbour][KEY] < g_bu[neighbour]):
            #     VLIST[neighbour][KEY] = g_bu[neighbour]
            #     heapq.heapify(VLIST)

    # base_ver = vertices[base_vertex_index]
    # create array of vertices in the area

    points_in_area_length = 0
    points_in_area = -np.ones(2000, dtype = np.int32)

    #Check if vertex is within area bounds
    for vertex_index, distance in enumerate(g_bu):
        if(distance <= r_threshold):
            #check if it has already been included in the area
            if (unvisited_vertices[vertex_index] and 
                vertex_index != base_vertex_index):
                #check the sizing
                if(points_in_area_length >= len(points_in_area)):
                    points_in_area = np.resize(points_in_area, (len(points_in_area) + 500,))
                    # npt.expand_base_mones(points_in_area, 500)
                points_in_area[points_in_area_length] = vertex_index
                points_in_area_length+= 1

            unvisited_vertices[vertex_index] = 0

    #Cant say Im sure as to why but in their implementation they repeat
    # this line here.        
    if(points_in_area_length >= len(points_in_area)):
        points_in_area = np.resize(points_in_area, (len(points_in_area) + 500,))
        # npt.expand_base_mones(points_in_area, 500)

    points_in_area[points_in_area_length] = base_vertex_index
    points_in_area_length+= 1
    # 0 is that it is no longer 'unvisited'
    unvisited_vertices[base_vertex_index] = 0 

   
    base_area_timer_start = time.process_time()
    calculateBaseArea(points_in_area, points_in_area_length, vertices,
                        triangles,  base_areas, base_points_length, A_matrix, debug_counters)
    base_area_timer_end = time.process_time()
    base_area_timer += base_area_timer_end - base_area_timer_start
    return g_bu, heap_restructuring_time, base_area_timer



def calculateBaseArea(points_in_area, points_in_area_length, vertices, triangles, 
                       base_areas, base_points_length, A_matrix, debug_counters):
    area = 0
    
    vertlen = len(vertices)
    for pointIndex in points_in_area:
        indexA = pointIndex
        neighbours_A = A_matrix[indexA, :].toarray()
        neighbours_A = neighbours_A.nonzero()[1]
        # neighbours_A = 
        for neighbour_A in neighbours_A:
            indexB = neighbour_A

            neighbours_B = A_matrix[indexB, :].toarray()
            neighbours_B = neighbours_B.nonzero()[1]
            for neighbour_B in neighbours_B:
                indexC = neighbour_B

                if (indexA != indexB and indexB != indexC and indexC != indexA):
                    area = area + calculateTrigArea(vertices[indexA], vertices[indexB], vertices[indexC])

    if (points_in_area_length == 0):
        area = 0
        debug_counters[0] += 1
    if (points_in_area_length == 1):
        area = 1
        debug_counters[1] += 1
        if(debug_counters[1] % 10 == 0):
            print(debug_counters[1])
        #At one point early in development this was a great showcase of logical errors.
    if (points_in_area_length == 2):
        area = 2
        debug_counters[2] += 1
    
    
    base_areas[base_points_length] = area



def calculateTrigArea(A, B, C):

    a = np.sqrt(np.power(C[0] - B[0], 2)
                +np.power(C[1] - B[1], 2)
                +np.power(C[2] - B[2], 2))
    
    b = np.sqrt(np.power(A[0] - C[0], 2)
                +np.power(A[1] - C[1], 2)
                +np.power(A[2] - C[2], 2))
    
    c = np.sqrt(np.power(A[0] - B[0], 2)
                +np.power(A[1] - B[1], 2)
                +np.power(A[2] - B[2], 2))
    
    p = (a + b + c) / 2
    area = (p *(p - a) * (p - b) * (p - c))

    if (area < 0):
        area = 1.0

    area = np.sqrt(area)

    return area


def decrease_key(VLIST, neighbour_index, distance_check):
    KEY = 0
    INDEX = 1
    vlist_index = 0
    for item in VLIST:
        if (item[INDEX] == neighbour_index):
            if(VLIST[vlist_index][KEY] > distance_check):
                VLIST[vlist_index][KEY] = distance_check
                heapq._siftup(VLIST, vlist_index)
                # heapq.heapify(VLIST)
            break
        vlist_index += 1
    return

def calculate_mu(g_values, base_points_length, point_index, base_areas):
    value = 0 
    for base_i in range(base_points_length):
        value = value + g_values[base_i][point_index] * base_areas[base_i]
    return value

def normalize_mu():
    pass


if __name__ == "__main__":

    create_test_mesh()