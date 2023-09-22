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
    
    resample_inbetweens(vertices, A_matrix, mu_values, ranges)


def resample_inbetweens(vertices, A_matrix, mu_values, ranges):
    for vertex_index in range(len(vertices)):
        range1 = gather_range(vertex_index, mu_values, ranges)
        neighbours = A_matrix[vertex_index, :].toarray()
        neighbours = neighbours.nonzero()[1]
        for neighbour in neighbours:
            range2 = gather_range(neighbour, mu_values, ranges)
            if(range1 != range2):
                if(range1 > range2):
                    new_mu = ranges[range2 + 1]
                else:
                    new_mu = ranges[range1 + 1]
                
                create_new_vertex(vertices, vertex_index, neighbour, mu_values, new_mu)
    return

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
    
    #Then you need to resize the vertices, and change the triangles accordingly
    #Good Fucking Luck with that, buddy.
    


#nice little numpy handling of the original implementation.
#Gets the mu range of input vertex
def gather_range(point_index, mu_values, ranges):
    if np.isclose(mu_values[point_index], 1, rtol = 1e-8, atol = 1e-6):
        return len(ranges) - 2
    left_range = ranges[0:-1]
    right_range = ranges[1:]
    vertex_rangeL = left_range <= mu_values[point_index]
    vertex_rangeR = mu_values[point_index] < right_range 
    vertex_range = np.where(vertex_rangeR * vertex_rangeL == 1)[0]
    return vertex_range


    


    

if __name__ == "__main__":
    gather_range(8, np.arange(10), np.arange(10))