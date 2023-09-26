import open3d as o3d
import open3d.visualization.gui as gui
import open3d.visualization.rendering as rendering
from open3d.visualization.gui import MouseEvent, KeyEvent
from open3d.visualization.rendering import Camera
import time 
import kmapper as km
import sklearn
from scipy.spatial.transform import Rotation as R
import numpy as np
# import kd_tree
# import hdbscan
# import trimesh as tm
# from matplotlib import pyplot as plt
from scipy.optimize import linprog
import copy


red = np.array([1,0,0])
green = np.array([0,1,0])
blue = np.array([0,0,1])
yellow = np.array([1,1,0])
magenta = np.array([1,0,1])
cyan = np.array([0,1,1])
black = np.array([0.1,0.1,0.1])
white = np.array([0.8,0.8,0.8])

def mesh_segmentation(vertices, geometry, t_param , resolution , convexity_measure, save_graphs = False):


    color_list = [red, green, blue, yellow, magenta, cyan, black]


     #Comforming with paper parameters, t_param = 33 would achieve the paper's accuracy.
    unit_sphere = o3d.geometry.TriangleMesh.create_sphere(resolution = 120)
    unit_sphere_sampling = np.asarray(unit_sphere.sample_points_uniformly(t_param).points)
    morse_directions = np.array([[0, 1, 0],[1, 0, 0],[0, 0, 1]])
    morse_directions = np.concatenate((morse_directions, unit_sphere_sampling), axis = 0)

    morse_functions = np.dot(vertices, morse_directions.T).T

    reeb_graphs = []
    reeb_graph_adjacencies = []
    mp_set = []
    candidate_cut_set = []
    
    mapper = km.KeplerMapper(verbose = 0)
    graph_acc_timer = 0
    mpis_acc_timer = 0
    cuts_acc_timer = 0
    # convexity_measure = 80
    # convexity_measure_nodes = 18 #default 30 on Armadillo
    debug_counter = 0
    eps = 0.14
    min_samples = 40
    perc_overlap = 0.12
    for i, morse_function in enumerate(morse_functions):
        graphs_acc_iter = time.process_time()
        graph = mapper.map(morse_function, vertices, #default eps = 12 on armadillo
                                    clusterer=sklearn.cluster.DBSCAN(eps=eps, min_samples=min_samples),
                                      cover = km.Cover(n_cubes = resolution, perc_overlap = perc_overlap))
        graph_acc_timer += time.process_time() - graphs_acc_iter
        reeb_graphs.append(graph)
        mpi, cci, rb_adj_i = iteration_mutex_pairs(graph, convexity_measure, debug_counter)
        mp_set.append(mpi)
        candidate_cut_set.append(cci)
        reeb_graph_adjacencies.append(rb_adj_i)
        if(save_graphs):
            mapper.visualize(graph, path_html = 'graphs/armadilloTests' + str(i) + '.html', title = 'eps:' + str(eps) +' ms:'+ str(min_samples)+' po:'+ str(perc_overlap))
        print(f'Computed Graph and Mutex Pairs of {i+1}/{t_param + 3}')
    print(f'Total Processing Time for Reeb Graphs: {graph_acc_timer}')




    ##COMPUTE COSTS AND SATISFACTION MATRIX##

    fcut_timer_start = time.process_time()
    final_cut_set = compute_final_cut_set(reeb_graphs, reeb_graph_adjacencies, mp_set, candidate_cut_set)
    fcut_timer_end = time.process_time() - fcut_timer_start
    print(f'Cut Set Calc Time: {fcut_timer_end}')

    cuts_cost_timer = time.process_time()
    objective_function = cost_of_cuts(reeb_graphs, candidate_cut_set, morse_directions, geometry)
    print(f'Cuts Cost Time: {time.process_time() - cuts_cost_timer}')
    # mapper.visualize(graph_x, path_html = 'armadilloX.html', title = 'eps0.06,ms = 35')

    A_ub = -np.array(final_cut_set).transpose()

    rhs = -np.ones(len(A_ub))
    # print('ok!')

    x_bounds = []
    for i in range(len(objective_function)):
        x_bounds.append([0, 1])
    
    for iter, line in enumerate(A_ub):
        if sum(line) == 0:
            print(iter)
    

    res = linprog(objective_function, A_ub = A_ub, b_ub = rhs,
              bounds = x_bounds)
    print(res)



    #for every direction, morse function is
    #the dot product of all points of the mesh/pointcloud
    #with that direction.
    #THEN, for every morse function (33 in number)
    #You calculate Reeb graph, which now is 
    #An one dimensional entity(?)

    # ##COLOR THE MESH's CUTS########
    # selected_cuts = res.x
    # print(selected_cuts)
    # colors = np.ones((len(vertices), 3)) * 0.6
    # colored_clusters = 0
    # counter = 0
    # for iter, reeb_graph in enumerate(candidate_cut_set):
    #     for j_iter, cut in enumerate(reeb_graph):
            
    #         if (selected_cuts[counter] > 0.49):
    #             color_pick = color_list[int(np.random.choice(len(color_list), 1))]
    #             colors[reeb_graphs[iter]['nodes'][cut]] = color_pick
    #             colored_clusters += 1
    #         counter += 1
    # print(f'Now colored {colored_clusters} clusters!')


    print(color_list)
    # geometry.vertex_colors = o3d.utility.Vector3dVector(colors)

    return res,candidate_cut_set, color_list, reeb_graphs



def downsample(point_cloud, a):

    points = np.asarray(point_cloud.points)
    N = np.shape(points)[0]

    indices = np.arange(N)
    M = N // a
    indices = np.random.choice(indices, M, replace = False)

    points = points[indices,:]

    point_cloud.points = o3d.utility.Vector3dVector(points)
    return point_cloud, M

def get_center(point_cloud):

    points = np.asarray(point_cloud.points)
    
    center = np.sum(points,axis=0) / points.shape[0]

    return center

def unit_sphere_normalization(point_cloud):

    points = np.asarray(point_cloud.points)

    distances = np.sum(np.square(points),axis=1)

    points = points / np.sqrt(np.max(distances))

    point_cloud.points = o3d.utility.Vector3dVector(points)
    return point_cloud

def translate(point_cloud, translation_vec):

    points = np.asarray(point_cloud.points)

    points += translation_vec

    point_cloud.points = o3d.utility.Vector3dVector(points)

    return point_cloud

# Curvature Tricks

def pca_compute(data, sort=True):
    """
	SVD decomposition
    """
    average_data = np.mean(data, axis=0) 
    decentration_matrix = data - average_data
    H = np.dot(decentration_matrix.T, decentration_matrix)
    eigenvectors, eigenvalues, eigenvectors_T = np.linalg.svd(H)
    if sort:
        sort = eigenvalues.argsort()[::-1]
        eigenvalues = eigenvalues[sort]
    return eigenvalues

def caculate_surface_curvature(cloud, radius=0.003):
    points = np.asarray(cloud.points)
    kdtree = o3d.geometry.KDTreeFlann(cloud)
    num_points = len(cloud.points)
    curvature = []  
    for i in range(num_points):
        k, idx, _ = kdtree.search_radius_vector_3d(cloud.points[i], radius)
        neighbors = points[idx, :]
        w = pca_compute(neighbors)
        delt = np.divide(w[2], np.sum(w), out=np.zeros_like(w[2]), where=np.sum(w) != 0)
        curvature.append(delt)
    curvature = np.array(curvature, dtype=np.float64)
    return curvature




# Reeb Graph Stuff

def iteration_mutex_pairs(graph, convexity_measure, debug_counter):
    j = 0
    clusters = []
    for i in range(graph['meta_data']['n_cubes']):
        while("cube"+str(i)+"_"+"cluster"+str(j) in graph['nodes']):
            j += 1
        clusters.append(j)
        j = 0
    return mp_search(graph, clusters, convexity_measure, debug_counter)


def mp_search(graph, clusters, convexity_measure, debug_counter):

    # nodes = {node: False for node in graph['nodes']}
    # notFound = True
    candidate_cuts = []
    mutex_pairs = {}
    # adjacency_graph = np.eye(len(graph['nodes']))
    neighbours = {node: [] for node in graph['nodes']}
    for node1 in neighbours:
        for node2 in graph['links'].copy():
            if exists_edge(graph, node1, node2):
                #This feels dumb and a clever way should exist somewhere.
                if (node2 not in neighbours[node1]):
                    neighbours[node1].append(node2)
                if (node1 not in neighbours[node2]):
                    neighbours[node2].append(node1)
    
    for node in graph['nodes']:
        hcube = node[4: node.index('_')]
        if (clusters[int(hcube)] > 1):
            mpi = 1
            while (mpi < clusters[int(hcube)]): 
                start = node
                end = 'cube' + hcube + '_' + 'cluster' + str(mpi)
                if (start == end):
                    mpi += 1
                    continue
                # debug_counter += 1
                # if debug_counter == 100:
                #     print('critical junction')
                convexity_mpi, validated = validated_bfs_rev(neighbours, start, end)
                mpi += 1
                if((convexity_mpi > convexity_measure) and validated):
                    mutex_pairs[start +' - ' + end ] = (graph['nodes'][start][0], graph['nodes'][end][0])
                #mp_reps.append((graph['nodes'][start][0], graph['nodes'][end][0])) 
        if(len(neighbours[node]) > 2 ):
            # candidate_cuts.append(node)
            for neighbour in neighbours[node]:
                candidate_cuts.append(neighbour)
    return mutex_pairs, candidate_cuts, neighbours




def exists_edge(graph, node1, node2):
    if ((node2 in graph['links'][node1]) or (node1 in graph['links'][node2])):
        return True
    return False


def bfs(graph, start, end):
    #Time save if you search for all X in hypercube
    # at the same bfs pass.

    visited = set()
    frontier = []
    visited.add(start)
    frontier.append(start)
    for degree in range(1, len(graph)):
        nextFrontier = []
        for item in frontier:
            for neighbour in graph[item]:
                if neighbour not in visited:
                    visited.add(neighbour)
                    nextFrontier.append(neighbour)
                    if neighbour == end:
                        return (degree + 1)    
        frontier = nextFrontier
    return False

def validated_bfs_rev(graph, start, end):
    #Time save if you search for all X in hypercube
    # at the same bfs pass.

    visited = set()
    frontier = []
    visited.add(start)
    frontier.append(start)
    validity_acc = 0
    validated = False
    for degree in range(1, len(graph)):
        nextFrontier = []

        for item in frontier:
            validity_acc += 1
            # if (len(frontier) > 1):
            #     validated = True
            for neighbour in graph[item]:
                if neighbour not in visited:
                    visited.add(neighbour)
                    nextFrontier.append(neighbour)
                    if neighbour == end:
                        # return degree + 1, validated
                        if(degree  >= validity_acc):
                            return (degree + 1), False
                        else:
                            return (degree + 1), True
                            
        frontier = nextFrontier
    return False, False

def validated_bfs(graph, start, end):
    
    visited = set()
    frontier = []
    visited.add(start)
    frontier.append(start)
    validator_flag = []
    validated = False
    for degree in range(1, len(graph)):
        nextFrontier = []
       
        for item in frontier:
            
            for neighbour in graph[item]:
                if neighbour not in visited:
                    visited.add(neighbour)
                    nextFrontier.append(neighbour)
                    if (len(graph[neighbour]) > 2):
                        validator_flag.append(len(nextFrontier))
                        # validated = True
                    elif(len(validator_flag) != 0):
                        frontier_item_index = frontier.index(item)
                        if (frontier_item_index in validator_flag):
                            validator_flag[validator_flag.index(frontier_item_index)] = len(nextFrontier)
                    if (neighbour == end):
                        # return degree + 1, validated
                        if (nextFrontier.index(neighbour) in validator_flag):
                            validated = True
                        return degree + 1, validated
                        # else: 
                        #     return degree , validated  
        if(len(nextFrontier) == len(validator_flag)):
            validated = True    
        frontier = nextFrontier
    return False, False



def compute_final_cut_set(reeb_graphs, reeb_graph_adjacencies, mutex_pair_sets, candidate_cuts):
    final_cut_set = []
    for reeb_graph_iteration, cuts_set in enumerate(candidate_cuts):
        for cut in cuts_set:
            cut_index = check_cut_validity(reeb_graphs[reeb_graph_iteration],
                                            reeb_graph_adjacencies[reeb_graph_iteration],
                                              mutex_pair_sets, cut)
            final_cut_set.append(cut_index)
    return final_cut_set

def check_cut_validity(reeb_graph_of_cut, reeb_graph_neighbours, mutex_pair_sets, cut):
    '''Should return ones for the mps the cut satisfies.'''
    mpi_row_check = []
    for reeb_graph_iteration in mutex_pair_sets:
        for mutex_pair in reeb_graph_iteration.items():
            start, end = linear_binary_search(reeb_graph_of_cut, mutex_pair[1])
            if (start == False or end == False):
                mpi_row_check.append(0)
                continue
            mpi_row_check.append(bfs_check(reeb_graph_neighbours, start, end, cut))
    return mpi_row_check
            
    

def bfs_check(graph, start, end, cut):
        #Graph here means neighbours graph, *Ahem*.

    if (start == False or end == False):
        return 0
    if (start  == end ):
        return 0
    
    visited = set()
    frontier = []
    cut_flag = []
    visited.add(start)
    frontier.append(start)
    for degree in range(1, len(graph)):
        nextFrontier = []
        
        for item in frontier:
            
            for neighbour in graph[item]:
                if neighbour not in visited:
                    visited.add(neighbour)
                    nextFrontier.append(neighbour)
                    if (item == cut):
                        cut_flag.append(len(nextFrontier))
                    elif(len(cut_flag) != 0):
                        frontier_item_index = frontier.index(item)
                        if (frontier_item_index in cut_flag):
                            cut_flag[cut_flag.index(frontier_item_index)] = len(nextFrontier)
                    if (neighbour == end):
                        if (nextFrontier.index(neighbour) in cut_flag):
                            return 1
                        else: return 0
        if(len(nextFrontier) == len(cut_flag)):
            #Small Optimisation
            # If whole Frontier is cut's neighbours, then surely cut satisfies condition.
            return 1      
        frontier = nextFrontier
    return 0

def linear_binary_search(reeb_graph, targets):
    '''Checks the clusters iteratively and performs 
    binary search to find the cluster in which the element belongs.'''
    start_found = False
    end_found = False
    #This corrects clusters overlap bug.
    start, end = False, False
    for node in reeb_graph['nodes']:
        if (binary_search(reeb_graph['nodes'][node], targets[0])):
                start_found = True
                start = node
        if (binary_search(reeb_graph['nodes'][node], targets[1])):
                end_found = True
                end =  node
        if (start_found and end_found):
            break
    return start, end  


def binary_search(arr, x):
    low = 0
    high = len(arr) - 1
    mid = 0
 
    while low <= high:
 
        mid = (high + low) // 2
 
        # If x is greater, ignore left half
        if arr[mid] < x:
            low = mid + 1
 
        # If x is smaller, ignore right half
        elif arr[mid] > x:
            high = mid - 1
        # means x is present at mid
        else:
            return True

    return False



#Computational Geometry And Cost of Cuts.

def cost_of_cuts(reeb_graphs, candidate_cut_set, morse_directions, mesh):
    costs_of_cuts = []
    for reeb_graph, cuts in enumerate(candidate_cut_set):
        for cut in cuts:
            cost = compute_cost(reeb_graphs[reeb_graph]['nodes'][cut], morse_directions[reeb_graph], mesh)
            costs_of_cuts.append(cost)
    return costs_of_cuts


def compute_cost(cut, morse_direction, mesh):

    rotated_mesh = rotate_mesh(mesh, morse_direction)
    vertices = np.asarray(rotated_mesh.vertices)
    triangles = np.asarray(rotated_mesh.triangles)

    # plane_height = -np.dot(vertices[cut[0]], morse_direction) 
    plane_height = -vertices[cut[0], 1]

    # print(f'PH:{plane_height}')
    plane_vec = np.array([0, 1, 0, plane_height]) #get Y cord of element in cut cluster
    # print(plane_vec)
    intersect_idxs = find_mesh_plane_intersection(vertices, triangles, plane_vec)
    intersect_triangles = triangles[intersect_idxs, :]
    # vertices_preselected = vertices[intersect_triangles[:, 0]] #Take only one vertex of each triangle, its enough
    # print(f'len before: {len(intersect_triangles)}')

    cluster_vertices = []
    for vertex in intersect_triangles[:, 0]:
        if (binary_search(cut, vertex)):
            cluster_vertices.append(vertex)
    #from selected you need to use only those 
    #That are on the cut-cluster.
    # print(f'len after: {len(cluster_vertices)}')
    vertices_selected = vertices[cluster_vertices, :]
    if (len(vertices_selected) != 0):
        plane_projection = plane_vec[:3]
        dot_n = np.dot(vertices_selected, plane_projection)[:, np.newaxis] * plane_projection

        projection_2d = np.subtract(vertices_selected, dot_n)
        sorted_projection = pcd_anglesort(projection_2d)
        # polygon_area = vertices_to_lineset(sorted_projection)

        area = area_shoelace(sorted_projection)
        # print(f'polygon area: {area} squared units')
        return area
    # print('nothing found, nothing gained.')
    return 100000000 #intersection of cut not in cluster. So cost is Veeeery large number.
                    #If picked its rip.



def find_mesh_plane_intersection(vertices, triangles, plane_vec):
    # print(plane_vec)
    vertex_map = np.concatenate([vertices,\
                                  np.ones((vertices.shape[0], 1))], axis = -1)
    
    # get the number of triangles in the mesh
    num_triangles = triangles.shape[0]
    # unroll the vertices, so that continous triad of vertices 
    # will belong in the same triangle
    # calculating the distance of the points to the plane
    unrolled_distances_mapped = vertex_map.dot(plane_vec)
    unrolled_distances_mapped = unrolled_distances_mapped > 0
    # see if the distance is positive or negative
    # reshaping the triangles back to original shape
    #same as before, open mapped vertices so as to have the checks
    #in triangle-list compatible form. 
    open_triangles = unrolled_distances_mapped[triangles]
    open_triangles.reshape(-1, 3)
    # create a boolean array to store the indexes of the intersection 
    # initially all points are intersection candidates
    intersect_idxs = np.ones((num_triangles, ), dtype=bool)
    intersect_idxs[open_triangles.all(axis=-1)] = 0
    #triangles_lower = intersect_idxs.copy()
    open_triangles = np.logical_not(open_triangles)
    # triangles that all vertices are below the plane
    intersect_idxs[open_triangles.all(axis=-1)] = 0
    return intersect_idxs #, triangles_lower, triangles_upper




def pcd_anglesort(pointcloud):
    p0 = pointcloud[np.argmin(pointcloud[:,2])]

    angles = np.arctan2(pointcloud[:,2] - p0[2], pointcloud[:,0] - p0[0])

    # Sort the points by angle
    sorted_idx = np.argsort(angles)
    pointcloud = pointcloud[sorted_idx]
    return pointcloud

def vertices_to_lineset(points, color=black):
    #used for debugging projection of vertices onto 2d plane

    # points = pad_2d(points, None, -1, 0)

    indices1 = np.arange(points.shape[0])
    indices2 = np.arange(points.shape[0])+1
    indices2[-1] = 0

    indices = np.vstack((indices1, indices2)).T

    return o3d.geometry.LineSet(
        o3d.utility.Vector3dVector(points),
        o3d.utility.Vector2iVector(indices)
    ).paint_uniform_color(color)


#method for calculating area for 2d polygon
def area_shoelace(projections_3d):
    '''Input is 3d [xcords, random, zcords] sorted by xz angle.
    Returns area of irregular polygon.'''

    iter = np.arange(len(projections_3d))
    x_i = np.array(projections_3d[:, 0])
    y_i = np.array(projections_3d[:, 2])
    area = np.sum(x_i[iter-1]*y_i[iter]
                  -x_i[iter]*y_i[iter-1]) * 0.5
    return area

    
def rotate_mesh(mesh, direction):
    rotated_mesh = copy.deepcopy(mesh)
    axis = np.cross(direction, [0, 1, 0])
    angle = np.arccos(np.dot(direction, [0, 1, 0]))
    rot_vector = angle * axis
    R = o3d.geometry.get_rotation_matrix_from_axis_angle(rot_vector)
    rotated_mesh.rotate(R)
    return rotated_mesh





