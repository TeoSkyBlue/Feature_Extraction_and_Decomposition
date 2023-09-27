from rb_graph_class import *


def compare_reeb_graphs(mu_coeff, MRG_num, points_num, weight, graph1, graph2):
    SIM_R_S = 0
    attributes1 = None
    attributes2 = None
    MRG1 = graph1
    MRG2 = graph2

    #Maybe the MRG's should be read from file instead of being generated
    #at runtime
    if (len(MRG1)!= len(MRG2)):
        print("The graphs have different resolution, aborting comparison.")
        return
    
    #calculate the attributes for each graph
    assignRestAttribs(MRG1, attributes1)
    assignRestAttribs(MRG2, attributes2)
    
    #compute parents in each graph MRG
    MRG1, attributes1, MRG2, attributes2 = computeParents(MRG1, attributes1, MRG2, attributes2)   

    SIM_R_S = make_comparison(MRG1, attributes1, MRG2, attributes2, mu_coeff, MRG_num, points_num, weight, SIM_R_S)

    return SIM_R_S



def make_comparison(MRG1, attributes1, MRG2, attributes2, mu_coeff, MRG_num, points_num, weight, SIM_R_S):
    NLIST = []
    MPAIR = []
    list1 = []
    list2 = []

    #Put nodes in from coarsest RES into NLIST
    reeb_graph = MRG1[-1]
    for i in range(len(reeb_graph)):
        vec = reeb_graph[i]
        node = Rnode(vec[0].index, vec[0].parent, vec[0].attribute.a, vec[0].attribute.l)
        list1.append(node)
    
    reeb_graph = MRG2[-1]
    for i in range(len(reeb_graph)):
        vec = reeb_graph[i]
        node = Rnode(vec[0].index, vec[0].parent, vec[0].attribute.a, vec[0].attribute.l)
        list2.append(node)

    NLIST.append(list1)
    NLIST.append(list2)
    
    while(len(list1) != 0 and len(list2) != 0):
        look_for_matching_pair(list1, list2, MPAIR, mu_coeff, MRG_num, points_num, weight, SIM_R_S)
    
    #total similarity 
    for i in range(len(MPAIR)):
        el = MPAIR[i]
        SIM_R_S = SIM_R_S + sim(el.node1, el.node2, weight)

    print("This passing a simple smoke test makes me happy.")
    print("More than you believe.")
    print("What do you mean, similarity?")
    print("\n\n\n..Oh, that too. It is " + str(SIM_R_S) + "between these two.")

    return SIM_R_S


def sim(node1, node2, weight):
    if(node1.attribute.a > node2.attribute.a):
        min_a = node2.attribute.a
    else:
        min_a = node1.attribute.a

    if(node1.attribute.l > node2.attribute.l):
        min_l = node2.attribute.l
    else:
        min_l = node1.attribute.l
    
    return (weight * min_a + (1 - weight) * min_l)


# takes node from NLIST with the largest sim(m, m) value and finds a pair for it
def look_for_matching_pair(NLIST, list1, list2, MPAIR, mu_coeff, MRG_num, points_num, weight, SIM_R_S):
    list1 = NLIST[0]
    list2 = NLIST[1]

    if(len(list1) == 0 or len(list2) == 0):
        return
    

    #Variable initialization

    is_first_time = True
    select_nlist2 = -1
    index2 = -1
    same_range_node_num = None
    select_nlist = None 
    index = None
    pair = MPairElement()
    tmp_v = findMaximumSim() #returns node with largest sim(m, m) value
    pair = tmp_v[0]
    select_nlist_o = tmp_v[1]
    index_o = tmp_v[2]
    nlist = select_nlist_o
    index = index_o

    n_list = None
    vec_m = None
    vec_n = None
    node = None
    matched_node = None
    temp_node = None
    m = None
    n = None

    #A bunch of other temps but python doesnt do tempts, its for .
    what_res_m = -1
    what_res_n = -1

    max_mat = -100000000000

    #Now look for a matching pair
    if(select_nlist == 2):
        n_list = NLIST[0]
        node = pair.node2
    else:
        n_list = NLIST[1]
        node = pair.node1
    first = True
    check = False
    for i in range(len(n_list)):
        temp_node = n_list[i]
        if(temp_node.left_bound == node.left_bound and
           temp_node.right_bound == node.right_bound):
            if (select_nlist == 2):
                m = temp_node
                n = node
            else:
                m = node
                n = temp_node

            range_diff = m.right_bound - m.left_bound
            what_res_m = calculateIndexInMRG1(MRG1, m)
            what_res_n = calculateIndexInMRG2(MRG2, n)
            



def computeParents(MRG1, attributes1, MRG2, attributes2):
    element, element2 = Rnode(), Rnode()
    for i in range(len(MRG1)):
        reeb_graph = MRG1[i]
        for j in range(len(reeb_graph)):
            vec = reeb_graph[j]
            element = vec[0]

            if (i == len(MRG1)-1):
                element.parent = None
            if(element.children != None):
                for k in range(len(element.children)):
                    temp_int = element.children[k]
                    next_reeb_graph = MRG1[i - 1]
                    vec2 = next_reeb_graph[temp_int]
                    element2 = vec2[0]
                    element2.parent = element

    for i in range(len(MRG2)):
        reeb_graph = MRG2[i]
        for j in range(len(reeb_graph)):
            vec = reeb_graph[j]
            element = vec[0]

            if(element.children != None):
                for k in range(len(element.children)):
                    temp_int = element.children[k]
                    next_reeb_graph = MRG2[i - 1]
                    vec2 = next_reeb_graph[temp_int]
                    element2 = vec2[0]
                    element2.parent = element

    # The way these returns work is unclear as to the way the objects are 
    #eventually transformed to Rnodes instead of ReebGraphs.
    return MRG1, attributes1, MRG2, attributes2


def assignRestAttribs(MRG, attributes):
    reeb_graph = MRG[0]
    for i in range(len(reeb_graph)):
        elements = reeb_graph[i]
        Rnode = Rnode(elements[0].index, elements[0].parent, attributes[i].a, attributes[i].l)
        #The way these objects interact looks fishy at best.
    
    for i in range(1, len(MRG)):
        reeb_graph = MRG[i]
        for j in range(len(reeb_graph)):
            elements = reeb_graph[j]
            Rnode = Rnode(elements[0].index, elements[0].parent, attributes[i].a, attributes[i].l)
            prev_graph = MRG[i-1]

            for k in range(len(Rnode.children)):
                temp_int = Rnode.children[k]
                adj_vec = prev_graph[temp_int]
                temp_node = adj_vec[0]
                Rnode.attribute.a = Rnode.attribute.a + temp_node.attribute.a
                Rnode.attribute.l = Rnode.attribute.l + temp_node.attribute.l
    
    return MRG, attributes

