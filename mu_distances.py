from utility import *


def geodesic_dijkstra(vertices, triangles, S_area):
    
    INF = 100000
    vertlen = len(vertices)
    r_threshold = np.sqrt(0.01 * S_area)
    print("r:", r_threshold)
    A_matrix = adjacency_matrix(triangles)
    # base_areas = np.zeros(10)
   
    unvisited_vertices = np.ones(vertlen, dtype = np.int8) #idxs of visited vertices
    base_default = 200

    g_values = -np.ones((base_default, vertlen), dtype = np.float128)
    mu_values = np.empty(vertlen, dtype = np.float128)
    
    

    base_points = -np.ones(base_default, dtype = np.int32)
    base_points_length = 0
    base_areas = -np.ones(base_default, dtype = np.float64)
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
    calculateBaseArea_co(points_in_area, points_in_area_length, vertices,
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


def calculateBaseArea_co(points_in_area, points_in_area_length, vertices, triangles, 
                       base_areas, base_points_length, A_matrix, debug_counters):
    area = 0
    

    # flat_triangles = triangles.flatten()
    # triangles_in_area =  np.intersect1d(flat_triangles, np.where(points_in_area >= 0))
    valid_points_in_area = points_in_area[np.where(points_in_area >= 0)]
    # mask = np.all(np.isin(triangles, valid_points_in_area), axis=0)
    #I want to find in which triangles do the points in the area exist:
    mask = np.isin(triangles, valid_points_in_area)

    triangles_in_area_mask = np.any(mask, axis = 1)
    
    triangles_in_area = triangles[triangles_in_area_mask]
    # vertices_in_area = vertices[valid_points_in_area]
    # triangles_in_area_flat = triangles_in_area.flatten()
    # vertices_in_area = vertices[np.unique(triangles_in_area_flat)]
    area = calculateAreaFromTriangles(triangles_in_area, vertices)


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


def calculateAreaFromTriangles(triangles_in_area, vertices_in_area):
    area_geometry = o3d.geometry.TriangleMesh(o3d.utility.Vector3dVector(vertices_in_area), o3d.utility.Vector3iVector(triangles_in_area))
    area = area_geometry.get_surface_area()
    return area
#o3d.geometry.TriangleMesh(o3d.utility.Vector3dVector(vertices_in_area), o3d.utility.Vector3iVector(triangles_in_area)).get_surface_area()

    



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