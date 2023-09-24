from mu_distances import *
from rb_graph_class import *

#Construct a multiresolutional Reeb Graph given the mesh and its 
#mu values
def geodesic_reeb_graph_extraction(mu_values, A_matrix, vertices, triangles, mrg_num):
    FINEST_RES = mrg_num
    percentage_threshold = 0.01
    all_t_sets = []
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
    # This kinda sucks, I know, sorry.
    triangles_to_keep_mask, new_triangles, new_vertices, new_triangles_counter, new_vertices_counter = resample_inbetweens(vertices, triangles, A_matrix, mu_values, ranges)
    #add new vertices to data structure
    vertices_rsd = np.vstack((vertices, new_vertices[0:new_vertices_counter]))
    #remove obsolete triangles from data structure
    triangles = triangles[triangles_to_keep_mask]
    #add the new triangles
    triangles_rsd = np.vstack((triangles, new_triangles[0:new_triangles_counter]))

    print("Number of triangles changed from ",old_triangle_num, " to ", len(triangles_rsd))
    print("Number of vertices changed from ",old_vertices_num, " to ", len(vertices_rsd))

    A_matrix_rsd = adjacency_matrix(triangles)

    #Missing checker logic to ensure intended behaviour in 
    #resampling.

    create_Tsets()



def create_Tsets(vertices, mu_values, ranges, FINEST_RES):
    point_ranges = -np.ones(len(vertices), dtype = np.float32)
    info = np.ones(len(vertices), dtype = np.int8)
    edged_counter = 0
    for i in range(len(vertices)):
        range_slot = gather_range(i, mu_values, ranges)
        range = ranges[range_slot]
        mu_value = mu_values[i]
        if(range == mu_value and range != 0):
            info[i] = 2
            edged_counter += 1

    lens = np.zeros(FINEST_RES)
    i = 0 
    while (i < len(info))

    #for all points
    #get point range
    #get mu_value
    #if range == muValue and range != 0
    #-> the point's info = 2
    
    #otherwise, its info is one


def resample_inbetweens(vertices, triangles, A_matrix, mu_values, ranges):
    

    triangles_to_keep_mask = np.ones(len(triangles), dtype = bool)
    # The 6 should be derived experimentally as to be optimal space complexity, 
    # no time for that at this moment.
    new_vertices = np.ones((6 *len(vertices),3), dtype = np.float32)
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
                    new_muAB = ranges[rangeB + 1]
                else:
                    new_muAB = ranges[rangeA + 1]
                vertexAB = create_new_vertex(vertices, indexA, indexB, mu_values, new_muAB)
                new_vertices[new_vertices_counter] = vertexAB
                vertex_indexAB = new_vertices_counter
                new_vertices_counter += 1
                
            if (checkerAC):
                if (rangeA > rangeC):
                    new_muAC = ranges[rangeC + 1]
                else:
                    new_muAC = ranges[rangeA + 1]
                vertexAC = create_new_vertex(vertices, indexA, indexC, mu_values, new_muAC)
                new_vertices[new_vertices_counter] = vertexAC
                vertex_indexAC = new_vertices_counter
                new_vertices_counter += 1
                
            if (checkerCB):
                if (rangeB > rangeC):
                    new_muBC = ranges[rangeC + 1]
                else:
                    new_muBC = ranges[rangeB + 1]
                vertexBC = create_new_vertex(vertices, indexB, indexC, mu_values, new_muBC)
                new_vertices[new_vertices_counter] = vertexBC
                vertex_indexBC = new_vertices_counter
                new_vertices_counter += 1

            new_triangles[new_triangles_counter : new_triangles_counter + 3] = alter_connections(vertex_indexAB,
                               vertex_indexBC, vertex_indexAC, triangle, checkerAB, checkerAC)
            new_triangles_counter += 3
            
    return triangles_to_keep_mask, new_triangles, new_vertices, new_triangles_counter, new_vertices_counter

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
    return np.array((p_x[0], p_y[0], p_z[0]))
    

    


#nice little numpy handling of the original implementation.
#Gets the mu range of input vertex
def gather_range(point_index, mu_values, ranges):
    #relative tolerance should be way lower than atol 
    #when handling values that go near 0
    if np.isclose(mu_values[point_index], 1, rtol = 1e-8, atol = 1e-6):
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