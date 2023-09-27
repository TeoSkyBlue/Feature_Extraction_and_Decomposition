from mu_distances import *
from rb_graph_class import *
# import kmapper as km
# import sklearn
from rb_graph_class import *
from copy import deepcopy


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

    # create_Tsets()
    # createFinestResReebGraph()
    # createMRG()
    #calculateAttributes()


    return #a whole lot of stuffs mainly: attributes, MRG, reebs, triangles_rsd, vertices_rsd, all_Tsets, A_matrix_rsd



def calculateAttributes(vertices, triangles, mu_values, ranges, FINEST_RES, all_Tsets, A_matrix, reebs, MRG, S_area):
    attributes = [attributeElement() * len(reebs[0])]
    lens = [None * len(all_Tsets)]
    rnum = len(MRG)
    sum_len = 0
    length = 0
    for i in range(len(all_Tsets)):
        Tset = all_Tsets[i]
        temp = attributeElement()
        #compute value of a(m) for each Tset
        temp.a = (1/rnum) * (calculateTsetArea(Tset, vertices, triangles) / S_area)
        attributes[i] = temp
        min_m = np.min(mu_values[Tset])
        max_m = np.max(mu_values[Tset])
        length = max_m - min_m
        sum_len += length
        lens[i] = length
    
    #Calculate l(m) of each Tset
    for i in range(len(all_Tsets)):
        attributes[i].l = (1 / rnum) * (lens[i] / sum_len)


    return attributes

#given that Tset is just a list of triangle indexes, we can calculate Tset area the same
# way we calculated the base areas on mu_approx.
def calculateTsetArea(Tset, vertices, triangles):
    Tset_area_geometry = o3d.geometry.TriangleMesh(o3d.utility.Vector3dVector(vertices), o3d.utility.Vector3iVector(triangles[Tset]))
    Tset_area = Tset_area_geometry.get_surface_area()
    return Tset_area




def createMRG(vertices, mu_values, ranges, FINEST_RES, all_Tsets, A_matrix, reebs, MRG):
    # main loop that creates all other resolutions
    prev_MRG_count = 0
    current_MRG_count = 0

    while(len(ranges) > 2):
        current_MRG_count = prev_MRG_count + 1
        #You are hoping for this to be copying a whole row of reeb class objects
        #Bold move Cotton, lets see if it pays out for him.
        reebs[current_MRG_count] = deepcopy(reebs[prev_MRG_count])

        #Start unifying the reeb graphs
        #Cons holds info in nodes to unify
        cons = [len(reebs[prev_MRG_count]) * None]
        for i in range(len(cons)):
            vec = MRG[prev_MRG_count][i]
            cons[i] = np.zeros(len(vec))

            for j in range(vec[0]):
                cons[i][j] = vec[j]

        for i in range(0, len(ranges) - 2, 2):
            unifyTwoRanges(i, i + 1, i + 2, prev_MRG_count, ranges, reebs, cons)
        
        #remove all other elements from ranges
        temp_ranges = ranges
        ranges = [None * (int(len(ranges) / 2) + 1)]
        j = 0
        for i in range(0, len(ranges), 2):
            ranges[j] = temp_ranges[i]
            j += 1
        ranges[0] = temp_ranges[0]
        ranges[len(ranges) - 1] = temp_ranges[len(temp_ranges) - 1]

        updateReebs(prev_MRG_count, reebs, MRG, cons)
        
        updateMRG(prev_MRG_count, reebs, MRG, cons)

        counter = 0
        for i in range(len(MRG[current_MRG_count])):
            el1 = reebs[current_MRG_count][i]
            adj = MRG[current_MRG_count][i]
            for j in range(1, adj[0]):
                intex = adj[j]
                el2 = reebs[current_MRG_count][intex]
                if(el1.left_bound == el2.left_bound and el1.right_bound == el2.right_bound):
                    counter += 1
        prev_MRG_count += 1
    
    return MRG, reebs

def updateReebs(prev_MRG_index, reebs):
    curr_MRG_index = prev_MRG_index + 1
    null_count = 0
    temp_reeb = reebs[curr_MRG_index]

    for i in range(len(reebs[curr_MRG_index])):
        if(temp_reeb == None):
            null_count += 1
    
    reebs[curr_MRG_index] = [None * (len(reebs[prev_MRG_index]) - null_count)]
    null_count = 0
    for i in range(len(temp_reeb)):
        if(temp_reeb[i] != None):
            temp_reeb[i].index = null_count
            reebs[curr_MRG_index][null_count] = temp_reeb[i]
            null_count += 1
    
    return reebs


def updateMRG(prev_MRG_index, reebs, MRG, cons):
    curr_MRG_index = prev_MRG_index + 1
    size_reeb = len(reebs[curr_MRG_index])
    parents = [None * size_reeb]

    for i in range(size_reeb):
        pars = reebs[curr_MRG_index][i].parents
        parents[i] = [len(pars) * None]
        
        for j in range(len(pars)):
            parents[i][j] = pars[j]

    
    MRG[curr_MRG_index] = [None * size_reeb]

    for i in range(size_reeb):
        MRG[curr_MRG_index][i].append(1)
    
    for i in range(size_reeb):
        j = i + 1
        while (j < size_reeb):
            for k1 in range(len(parents[i])):
                for k2 in range(len(parents[j])):

                    if(areNodesConnected(prev_MRG_index, parents[i][k1], parents[j][k2]) == True 
                       and areNodesConnected(curr_MRG_index, i, j) == False):
                        index1 = i
                        index2 = j
                        len1 = MRG[curr_MRG_index][index1][0]
                        len2 = MRG[curr_MRG_index][index2][0]

                        adj1 = MRG[curr_MRG_index][index1]
                        adj2 = MRG[curr_MRG_index][index2]
                        # Resizing adj1
                        # if(len1>= len(MRG[curr_MRG_index][index1])):
                        #     adj1 = np.hstack((adj1, np.zeros(len1)))
                        MRG[curr_MRG_index][index1][len1] = index2
                        MRG[curr_MRG_index][index1][0] += 1

                        # Resizing adj2
                        # if(len2>= len(MRG[curr_MRG_index][index2])):
                        #     adj2 = np.hstack((adj2, np.zeros(len2)))
                        MRG[curr_MRG_index][index2][len2] = index1
                        MRG[curr_MRG_index][index2][0] += 1
            j += 1
    return MRG


def areNodesConnected(mrg_index, index1, index2, MRG):
    adj = MRG[mrg_index][index1]
    len = adj[0]
    for i in range(1, len):
        if(adj[i] == index2):
            return True
    return False





def unifyTwoRanges(left_range_index, mid_range_index, right_range_index, prev_MRG_index, ranges, reebs, cons):
    element = ReebGraph()
    element2 = ReebGraph()

    curr_MRG_index = prev_MRG_index + 1
    left_range = ranges[left_range_index]
    mid_range = ranges[mid_range_index]
    right_range = ranges[right_range_index]

    i = 0 
    while ( i < len(reebs[curr_MRG_index])):
        if(reebs[curr_MRG_index][i] != None):
            element = reebs[curr_MRG_index][i]
            if((element.left_bound == left_range and element.right_bound == mid_range)
               or (element.left_bound == mid_range and element.right_bound == right_range)
               or (element.left_bound == left_range and element.right_bound == right_range)):
                j = 1 
                temp_i = i
                while (j < cons[element.index][0]):
                    index = cons[element.index][j]
                    if (reebs[curr_MRG_index][index] != None):
                        element2 = reebs[curr_MRG_index][index]
                        if((element2.left_bound == left_range and element2.right_bound == mid_range)
                            or (element2.left_bound == mid_range and element2.right_bound == right_range)
                            or (element2.left_bound == left_range and element2.right_bound == right_range)):
                            unifyTwoNodes(element.index, element2.index, left_range, right_range, prev_MRG_index, reebs, cons)
                            temp_i -= 1
                    j += 1      
                i = temp_i
        i += 1

    for i in range(len(reebs[curr_MRG_index])):
        if (reebs[curr_MRG_index][i] != None):
            elt = reebs[curr_MRG_index][i]

            if((elt.left_bound == left_range and elt.right_bound == mid_range)
            or (elt.left_bound == mid_range and elt.right_bound == right_range)):
                elt.left_bound = left_range
                elt.right_bound = right_range

                reebs[curr_MRG_index][i] = elt
    return #probably the structures you altered.
            

def unifyTwoNodes(node_index1, node_index2, left_b, right_b, prev_MRG_index, reebs, cons):
    curr_MRG_index = prev_MRG_index + 1
    bigger = 0
    smaller = 0
    if(node_index1 == node_index2):
        print("Warning, unifyTwoNodes reports that node_indexes 1 and 2 are the same.")
        return
    el1 = reebs[curr_MRG_index][node_index1]
    el2 = reebs[curr_MRG_index][node_index2]

    #here implementation uses size but that is just java's vector way of saying current len 
    #..I think
    for i in range(len(el2.Tsets)):
        temp_int = el2.Tsets[i]
        if(temp_int not in el1.Tsets):
            el1.Tsets.append(temp_int)

    for i in range(len(el2.parents)):
        int_t = el2.parents[i]
        if(int_t not in el1.parents):
            el1.parents.append(int_t)
    
    el1.left_bound = left_b
    el1.right_bound = right_b

    el2.left_bound = left_b
    el2.right_bound = right_b

    el2.index = -1

    reebs[curr_MRG_index][node_index2] = None
    adj1 = []
    len1 = 0
    adj2 = cons[node_index2]
    len2 = adj2[0]

    index1, index2 = 0, 0
    found = False

    for i in range(len2):
        index2 = adj2[i]
        adj1 = cons[node_index1]
        len1 = adj1[0]
        found = False

        for j in range(1, len1):
            index1 = adj1[j]
            if(index2 == index1):
                found == True
        
        if(found == False and index2 != node_index1):
            if(len1 >= len(adj1)):
                print("Resizing is needed on cons[node_index1]")
                pass
                # resize the cons[node_index1] array.
                cons[node_index1] 
            cons[node_index1][len1] = index2
            cons[node_index1][0] += 1
        
    return #the structures you altered.


def subrange_triangles(ranged_triangles, max_range, finest_res):
    ranges = np.arange(0, max_range)
    subranged_triangles = []
    low = 0
    for range_value in range(1, 1, len(ranges)):
        high = finest_res /(range_value + 1)
        subranged_triangles.append(np.where((ranged_triangles >= low) & (ranged_triangles <= high)))
        low = high
    return subranged_triangles

def create_Tsets(vertices, mu_values, ranges, FINEST_RES, all_Tsets, A_matrix):
    point_ranges = -np.ones(len(vertices), dtype = np.float32)
    info = np.ones(len(vertices), dtype = np.int8)
    border_counter = 0
    for i in range(len(vertices)):
        range_slot = gather_range(i, mu_values, ranges)
        rng = ranges[range_slot]
        mV = mu_values[i]
        if(rng == mV and rng != 0):
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
            one_set = create_one_Tset(i, rang, mu_values, ranges)
            one_range_Tsets = all_Tsets[rang]

            #resizing one_range_tsets

            one_range_Tsets[lens[rang]] = one_set
            lens[rang] += 1
            all_Tsets[rang] = one_range_Tsets


        elif (info[i] == 2 or info[i] == 3):
            rang = gather_range(i) - 1
            one_set = create_one_Tset(i, rang, mu_values, ranges)
            one_range_Tsets = all_Tsets[rang]

            #resizing one_range_tsets

            one_range_Tsets[lens[rang]] = one_set
            lens[rang] += 1
            all_Tsets[rang] = one_range_Tsets

    return all_Tsets #you can Trim the tsets if need be.


#     #for all points
#     #get point range
#     #get mu_value
#     #if range == muValue and range != 0
#     #-> the point's info = 2
    
#     #otherwise, its info is one

def create_one_Tset(point_index, range, vertices, info, mu_values, ranges, A_matrix):
    #will need to trim tset so that we dont have
    #gianormous data structures, but what can you do, really
    tset = np.zeros(len(vertices))
    tset_length = 0
    real_range = 0
    index = point_index
    stacker = []
    # temp = 0
    if(info[point_index] > 0):
        stacker.append(index)
    
    while(len(stacker) != 0):
        index = stacker.pop(0)
        #Resizing tset if needed
        stacker, point_result = create_one_point(index, range, mu_values, ranges, info, A_matrix, stacker)
        if(point_result == True):
        #     if(tset_length >= len(tset)):
        #         tset = np.hstack((tset, np.zeros(len(vertices))))
            #expansions of tset
            tset[tset_length] = index
            tset_length += 1

    #trim tset
    return tset


def create_one_point(point_index, range, mu_values, ranges, info, A_matrix, stacker):
    result = False
    real_range = gather_range(point_index, mu_values, ranges)
    if (info[point_index] > 0):
        result = True
        if(info[point_index] == 1):
            info[point_index] = 0
            if(real_range != range):
                print("Warning, somethings fishy in the ranges. 0")
        elif (info[point_index] == 2):
            if(real_range == range):
                info[point_index] = 3
            elif(real_range == range+1):
                info[point_index] = 1
            else:
                print("yea, that was not great either, warning on ranges. 1")
        
        elif(info[point_index] == 3):
            info[point_index] = 0
            if(real_range != range+1):
                print("warning, in ranges 2")

        neighbours = A_matrix[point_index, :].toarray()
        neighbours = neighbours.nonzero()[1]
        for neighbour in neighbours:
            if(info[neighbour] > 0):
                real_range = gather_range(neighbour, mu_values, ranges)
                if (info[neighbour] == 1 and real_range == range):
                    if(neighbour not in stacker):
                        stacker.append(neighbour)
                elif(info[neighbour] == 3 and real_range == range + 1):
                    if(neighbour not in stacker):
                        stacker.append(neighbour)
                elif(info[neighbour] == 2):
                    if(real_range == range + 1 or real_range == range):
                        if(neighbour not in stacker):
                            stacker.append(neighbour)
    return stacker, result
            

def createFinestResReebGraph(vertices, mu_values, ranges, A_matrix, triangles, all_Tsets):
    mrg_size = 0
    tm = len(ranges) - 1
    while(tm != 1):
        mrg_size += 1
        tm = tm / 2
    mrg_size += 1

    r_num = 0
    #calculate the number of nodes in the finest res reeb graph
    #all Tsets is 3 dimensional.
    for i in range(len(all_Tsets)):
        r_num += len(all_Tsets[i])
    
    MRG = [None * mrg_size] #int[MRG_size][][]
    MRG[0] = []
    MRG[0].append(np.ones(r_num, dtype = np.int32))

    reebs = [None * mrg_size]
    reebs[0].append([ReebGraph() * r_num])

    for i in range(r_num):
        rel = ReebGraph()
        rel.index = i
        temp_int = calculateRange(i, all_Tsets)
        rel.left_bound.append(ranges[temp_int])
        rel.right_bound.append(ranges[temp_int + 1])
        rel.Tsets.append(i)


        #Questionable logic, but it should work.
        reebs[0][0][i] = rel


    index1 = 0
    index2 = 0
    main_index = 0
    count = 0
    
    for i in range(len(all_Tsets) - 1):
        current = all_Tsets[i]
        following = all_Tsets[i + 1]

        for j in range(len(current)):
            Tset1 = current[j]
            for k in range(len(following)):
                Tset2 = following[k]
                if (isConnectedToTset(Tset1, Tset2, A_matrix)):
                    index2 = main_index + len(current) + k
                    index1 = main_index + j
                    len1 = MRG[0][index1][0]
                    len2 = MRG[0][index2][0]

                    adj1 = MRG[0][index1]
                    adj2 = MRG[0][index2]

                    #expansions of adj1
                    # if (len1 >= len(MRG[0][index1])):
                    #     adj1 = np.hstack((adj1, np.zeros(len1)))
                    MRG[0][index1][len1] = index2
                    MRG[0][index1][0] += 1


                    #expansions of adj2
                    # if(len2 >= len(MRG[0][index2])):
                    #     adj2 = np.hstack((adj2, np.zeros(len2)))

                    MRG[0][index2][len2] = index1
                    MRG[0][index2][0] += 1
        main_index += len(current)

    flatten_Tsets(all_Tsets)

    return MRG, reebs



def flatten_Tsets(all_Tsets):
    tcount = 0
    for i in range(len(all_Tsets)):
        tcount += len(all_Tsets[i])
    
    tsets = [None * tcount]
    k = 0
    for i in range(len(all_Tsets)):
        for j in range(len(all_Tsets[i])):
            tsets[k] = all_Tsets[i][j]
            k += 1

    return tsets



def isConnectedToTset(Tset1, Tset2, A_matrix):
    for i in range(len(Tset1)):
        temp = Tset1[i]
        if(isConnectedToTSet(temp, Tset2, A_matrix)):
            return True
    return False

def isConnectedToTSet(point_index, Tset, A_matrix):
    if(Tset == None or len(Tset) == 0):
        return True
    neighbours = A_matrix[point_index, :].toarray()
    neighbours = neighbours.nonzero()[1]

    for i in range(len(Tset)):
        index = Tset[i]
        for neighbour in neighbours:
            if (index == neighbour):
                return True
    
    return False





def calculateRange(point_index, all_Tsets):
    temp = 0
    for i in range(len(all_Tsets)):
        temp += len(all_Tsets[i])
        if(point_index < temp):
            return i
    
    print("Warning, calc range broky")
    return -1




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