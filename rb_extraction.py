from mu_distances import *
from rb_graph_class import *
# import kmapper as km
# import sklearn


#Construct a multiresolutional Reeb Graph given the mesh and its 
#mu values
def geodesic_reeb_graph_extraction(mu_values, A_matrix, vertices, triangles, base_vertices_num, r_threshold, mrg_num):
    FINEST_RES = mrg_num
    percentage_threshold = 0.01
    
    tsets = []
    
    ranges = []
    multires_rb_graph = dict()
    reeb_graphs = []
    r_num = 0
    
    #range division
    ranges = np.arange(0, FINEST_RES + 1)/ FINEST_RES
    ranges[FINEST_RES] = 1.0
    

    #Resampling Logic
    old_triangle_num = len(triangles)
    old_vertices_num = len(vertices)
    #Resampling Function
    # This kinda sucks, I know, sorry.
    triangles_to_keep_mask, new_triangles, new_vertices, new_triangles_counter, new_vertices_counter, mu_values = resample_inbetweens(vertices, triangles, A_matrix, mu_values, ranges)
    #add new vertices to data structure
    vertices_rsd = np.vstack((vertices, new_vertices[0:new_vertices_counter]))
    #remove obsolete triangles from data structure
    triangles = triangles[triangles_to_keep_mask]
    #add the new triangles
    triangles_rsd = np.vstack((triangles, new_triangles[0:new_triangles_counter]))

    print("Number of triangles changed from ",old_triangle_num, " to ", len(triangles_rsd))
    print("Number of vertices changed from ",old_vertices_num, " to ", len(vertices_rsd))
    # morse_function = mu_values

    A_matrix_rsd = adjacency_matrix(triangles)

    geometry_rsd = o3d.geometry.TriangleMesh(o3d.utility.Vector3dVector(vertices_rsd), o3d.utility.Vector3iVector(triangles_rsd))
    
    # for triangle in triangles_rsd:
    #     triangles_range_values = gather_range(mu_values, triangle[0])


    #A measly attempt to cheat the system based on past implementation, does not work
    # mapper = km.KeplerMapper(verbose = 0)
    # convexity_measure = FINEST_RES
    # graph = mapper.map(morse_function, vertices, #default eps = 12 on armadillo
    #                                 clusterer=sklearn.cluster.DBSCAN(eps=r_threshold, min_samples=1),
    #                                   cover = km.Cover(n_cubes = convexity_measure, perc_overlap = 0.2))
    # mapper.visualize(graph, path_html = 'armadilloTests.html', title = 'eps0.06,ms = 35')


    ranged_triangles = []
    triangles_in_range = np.empty(len(triangles_rsd))
    #range raming scheme is a tiny bit unfortunate, yes.
    for triangle in range(len(triangles_in_range)):
        triangles_in_range[triangle] = gather_range(triangles_rsd[triangle][0], mu_values, ranges)
    
    tset_geometries = []
    for range_value in range(len(ranges) - 1):
        ranged_triangles.append(np.where(triangles_in_range == range_value)[0])                                             #of the resampled triangles, get the ones that are on the range range_value
        tset_geometry = o3d.geometry.TriangleMesh(o3d.utility.Vector3dVector(vertices_rsd), o3d.utility.Vector3iVector(triangles_rsd[ranged_triangles[range_value]]))
        tset_geometries.append(tset_geometry)
    

    FRRG = [] 
    # subranged_triangles = subrange_triangles(triangles_in_range, 2, 8)
    # i = 1
    # while((2*i)<len(FINEST_RES)):
    #     subranged_triangles = subrange_triangles(ranged_triangles, i)
    #     temp_ranges = np.arange(1, 2*i)




    # Finest Res Reeb Graph (nodes only)
    for tset_range_group in tset_geometries:
        Rnodes = tset_range_group.cluster_connected_triangles()
        FRRG.append(Rnodes)
    

    #for each level,
    #   for earch node at that level
    #       chech every node above it and see if they are connected
    #if so, there exists an edge between those two.
    # FRRGE = []


    #This looks ugly. Probably because it is, but I can explain.
    #Ok, changed my mind, this is too ugly to explain.
    # for range_value in range(len(ranges) - 2):
    #     for tset in FRRG[range_value][0]:
    #         nodes_above = FRRG[range_value + 1][0]
    #         # concerned_triangles = triangles_rsd[nodes_above]
    #         for node_above in nodes_above:
    #             checker_triangles = np.vstack((triangles_rsd[node_above], triangles_rsd[tset]))
    #             checker_geometry = o3d.geometry.TriangleMesh(o3d.utility.Vector3dVector(vertices_rsd), o3d.utility.Vector3iVector(checker_triangles))
    #             _, clusters, _ = checker_geometry.cluster_connected_triangles()
    #             if (len(clusters) == 1):
    #                 FRRGE.append((tset, node_above))


    print("Finished will all Tsets")

    # all_Tsets = np.empty((FINEST_RES, 100, len(vertices_rsd)))
    #ALL_TSETS[FINEST_RES][number_of_sets][number_of_triangles_in_set]

    #Missing checker logic to ensure intended behaviour in 
    #resampling.

    #create_Tsets()



def subrange_triangles(ranged_triangles, max_range, finest_res):
    ranges = np.arange(0, max_range)
    subranged_triangles = []
    low = 0
    for range_value in range(1, 1, len(ranges)):
        high = finest_res /(range_value + 1)
        subranged_triangles.append(np.where((ranged_triangles >= low) & (ranged_triangles <= high)))
        low = high
    return subranged_triangles

def create_Tsets(vertices, mu_values, ranges, FINEST_RES, all_Tsets):
    point_ranges = -np.ones(len(vertices), dtype = np.float32)
    info = np.ones(len(vertices), dtype = np.int8)
    border_counter = 0
    for i in range(len(vertices)):
        range_slot = gather_range(i, mu_values, ranges)
        range = ranges[range_slot]
        mu_value = mu_values[i]
        if(range == mu_value and range != 0):
            info[i] = 2
            border_counter += 1
    print("Bordered vertices: ", border_counter)

    lens = np.zeros(FINEST_RES)
    i = 0 
    
    while (i < len(info)):
        if (info[i] == 0):
            i+= 1

        elif (info[i] == 1):
            rang = gather_range(i)
            one_set = create_one_Tset(i, rang)
            one_range_Tsets = all_Tsets[rang]

            #resizing one_range_tsets
            one_range_Tsets[lens[rang]] = one_set

            all_Tsets[rang] = one_range_Tsets


        elif (info[i] == 2 or info[i] == 3):
            rang = gather_range(i) - 1
            one_set = create_one_Tset(i, rang)
            one_range_Tsets = all_Tsets[rang]

            #resizing one_range_tsets
            one_range_Tsets[lens[rang]] = one_set

            all_Tsets[rang] = one_range_Tsets

    return all_Tsets #you can Trim the tsets if need be.


#     #for all points
#     #get point range
#     #get mu_value
#     #if range == muValue and range != 0
#     #-> the point's info = 2
    
#     #otherwise, its info is one

def create_one_Tset(point_index, range, vertices, info):
    #will need to trim tset so that we dont have
    #gianormous data structures, but what can you do, really
    tset = np.zeros(len(vertices))
    tset_length = 0
    real_range = 0
    index = point_index
    stacker = []
    if(info[point_index] > 0):
        stacker.append(index)
    
    while(len(stacker) != 0):
        temp = stacker.pop(0)
        if(create_one_point(index, range) == True):
            #expansions of tset
            tset[tset_length] = index
            tset_length += 1
    #trim tset
    return tset


def create_one_point():
    pass



def resample_inbetweens(vertices, triangles, A_matrix, mu_values, ranges):
    

    triangles_to_keep_mask = np.ones(len(triangles), dtype = bool)
    # The 6 should be derived experimentally as to be optimal space complexity, 
    # no time for that at this moment.
    new_vertices = np.ones((6 *len(vertices),3), dtype = np.float32)
    new_mu_values = np.zeros(6 *len(vertices), dtype = np.float32)

    new_triangles = np.ones((6 *len(triangles), 3), dtype = np.int32)
    new_triangles_counter = 0
    new_vertices_counter = 0
    for i, triangle in enumerate(triangles):
        indexA = triangle[0]
        rangeA = gather_range(indexA, mu_values, ranges)
        indexB = triangle[1]
        rangeB = gather_range(indexB, mu_values, ranges)
        indexC = triangle[2]
        rangeC = gather_range(indexC, mu_values, ranges)
        checkerAB = (rangeA != rangeB)
        checkerAC = (rangeA != rangeC)
        checkerCB = (rangeC != rangeB)
        vertexAB, vertexAC, vertexBC = 0, 0, 0
        vertex_indexAB, vertex_indexAC, vertex_indexBC,  = 0, 0, 0
        if (checkerCB or checkerAB or checkerAC):
            
            triangles_to_keep_mask[i] = 0
            if (checkerAB):
                if (rangeA > rangeB):
                    new_mu_values[new_vertices_counter] = ranges[rangeB + 1]
                else:
                    new_mu_values[new_vertices_counter] = ranges[rangeA + 1]
                vertexAB = create_new_vertex(vertices, indexA, indexB, mu_values, new_mu_values[new_vertices_counter])
                new_vertices[new_vertices_counter] = vertexAB
                vertex_indexAB = new_vertices_counter
                new_vertices_counter += 1
                
            if (checkerAC):
                if (rangeA > rangeC):
                    new_mu_values[new_vertices_counter] = ranges[rangeC + 1]
                else:
                    new_mu_values[new_vertices_counter] = ranges[rangeA + 1]
                vertexAC = create_new_vertex(vertices, indexA, indexC, mu_values, new_mu_values[new_vertices_counter])
                new_vertices[new_vertices_counter] = vertexAC
                vertex_indexAC = new_vertices_counter
                new_vertices_counter += 1
                
            if (checkerCB):
                if (rangeB > rangeC):
                    new_mu_values[new_vertices_counter] = ranges[rangeC + 1]
                else:
                    new_mu_values[new_vertices_counter] = ranges[rangeB + 1]
                vertexBC = create_new_vertex(vertices, indexB, indexC, mu_values, new_mu_values[new_vertices_counter])
                new_vertices[new_vertices_counter] = vertexBC
                vertex_indexBC = new_vertices_counter
                new_vertices_counter += 1

            new_triangles[new_triangles_counter : new_triangles_counter + 3] = alter_connections(vertex_indexAB,
                               vertex_indexBC, vertex_indexAC, triangle, checkerAB, checkerAC)
            new_triangles_counter += 3
    
    added_mu_values = np.hstack((mu_values, new_mu_values[0:new_vertices_counter]))
    return triangles_to_keep_mask, new_triangles, new_vertices, new_triangles_counter, new_vertices_counter, added_mu_values

def alter_connections(vertex_indexAB, vertex_indexBC, vertex_indexAC,
                       triangle, checkerAB, checkerAC):
    A, B, C = triangle[0], triangle[1], triangle[2]
    AB = vertex_indexAB
    BC = vertex_indexBC
    AC = vertex_indexAC

    if(checkerAB and checkerAC): #BC common range
        triangle1 = np.array([A, AB, AC])
        triangle2 = np.array([B, AB, C])
        triangle3 = np.array([C, AC, AB])
        
    elif(checkerAB and not checkerAC): #AC common range
        triangle1 = np.array([A, C, AB])
        triangle2 = np.array([C, BC, AB])
        triangle3 = np.array([B, AB, BC])
        
    elif(not checkerAB and checkerAC): #AB common range
        triangle1 = np.array([A, B, AC])
        triangle2 = np.array([B, AC, BC])
        triangle3 = np.array([C, AC, BC])

    return triangle1, triangle2, triangle3








def create_new_vertex(vertices, vertex_index, neighbour, mu_values, new_mu):
    point1 = vertices[vertex_index]
    point2 = vertices[neighbour]
    mu1 = mu_values[vertex_index]
    mu2 = mu_values[neighbour]
    if (mu1 > mu2):
        p_x = ((point1[0] * (new_mu - mu2)) + (point2[0] * (mu1 - new_mu))) / (mu1 - mu2)
        p_y = ((point1[1] * (new_mu - mu2)) + (point2[1] * (mu1 - new_mu))) / (mu1 - mu2)
        p_z = ((point1[2] * (new_mu - mu2)) + (point2[2] * (mu1 - new_mu))) / (mu1 - mu2)
    else:
        p_x = ((point1[0] * (mu2 - new_mu)) + (point2[0] * (new_mu - mu1))) / (mu2 - mu1)
        p_y = ((point1[1] * (mu2 - new_mu)) + (point2[1] * (new_mu - mu1))) / (mu2 - mu1)
        p_z = ((point1[2] * (mu2 - new_mu)) + (point2[2] * (new_mu - mu1))) / (mu2 - mu1)
    return np.array((p_x, p_y, p_z))
    

    


#nice little numpy handling of the original implementation.
#Gets the mu range index of input vertex
def gather_range(point_index, mu_values, ranges):
    #relative tolerance should be way lower than atol 
    #when handling values that go near 0
    if (np.isclose(mu_values[point_index], 1)):
        return len(ranges) - 2
    left_range = ranges[0:-1]
    right_range = ranges[1:]
    vertex_rangeL = left_range <= mu_values[point_index]
    vertex_rangeR = mu_values[point_index] < right_range 
    vertex_range = np.where(vertex_rangeR * vertex_rangeL == 1)[0]
    return vertex_range


    


    

if __name__ == "__main__":
    pass
    # gather_range(8, np.arange(10), np.arange(10))