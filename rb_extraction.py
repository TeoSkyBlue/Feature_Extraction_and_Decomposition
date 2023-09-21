from mu_distances import *

#Construct a multiresolutional Reeb Graph given the mesh and its 
#mu values
def geodesic_reeb_graph_extraction():
    FINEST_RES = 0
    percentage_threshold = 0.01
    all_t_sets = []
    tsets = []
    
    ranges = []
    multires_rb_graph = dict()
    reeb_graphs = []
    r_num = 0
    