- Feature_Extraction_and_Matching(not quite)
Part of coursework for 'Computational Geometry and Computer Vision' in 8th semester of Electrical and Computer Engineering dep. of University of Patras. 


 The project is now considered end of life, it is provided as is, should really not be used for anything that amounts to anything as it is largely incomplete.
 Within the repo exists some Feature Extraction functionality using python's open3d library
This includes: 
    - Pointcloud triangulation using existing open3d functionality
    - Convex hull
    - The visualization of curvature, delta coordinates, eigenvectors and 
    the integral of the geodesic distance
    - Isotropic Remeshing(forked from other repo by @sfcaracciolo)
    - Reeb Graph extraction
    - Nearest neighbour tree-based algorithms
    - Dijkstra approach on graph least distances
    - Mesh segmentation using the Reeb Graphs extracted, an approach 
    based on the work of Hairong Liu, Wenyu Liu, Longin Jan Latecki, 'Convex Shape Decomposition' (2010)

As you can see, there is no matching in the 'included' tab, as that was never finished, its files remain
in the project as a form of archaic artifacts of a time long gone.
Note that THERE ARE KNOWN ISSUES on the solutions provided, especially in the mesh segmentation implementation.

    Curvature visual on low poly armadillo

![alt text](https://github.com/TeoSkyBlue/Feature_Extraction_and_Matching/blob/main/screenshots/curvature.png?raw=true)	

    Delta Coordinates visual on low poly armadillo

![alt text](https://github.com/TeoSkyBlue/Feature_Extraction_and_Matching/blob/main/screenshots/delta.png?raw=true)

    Delta coordinates visual on armadillo after remeshing

![alt text](https://github.com/TeoSkyBlue/Feature_Extraction_and_Matching/blob/main/screenshots/deltaRemesh.png?raw=true)

    Integral of geodesic distance visual on low poly armadillo

![alt text](https://github.com/TeoSkyBlue/Feature_Extraction_and_Matching/blob/main/screenshots/geodesic.png?raw=true)

    Mesh segmentation results on original armadillo mesh/pointcloud 
    (method is vertex/point based, would work on either)

![alt text](https://github.com/TeoSkyBlue/Feature_Extraction_and_Matching/blob/main/screenshots/segmentation.png?raw=true)




    

    
