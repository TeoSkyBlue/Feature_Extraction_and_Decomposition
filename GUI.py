import open3d as o3d
import open3d.visualization.gui as gui
import open3d.visualization.rendering as rendering
from open3d.visualization.gui import MouseEvent, KeyEvent
from open3d.visualization.rendering import Camera
import sys
import os
import numpy as np
import utility as U
import curvature as crv
import platform
from scipy.sparse.linalg import eigs, eigsh
from scipy.linalg import eigh

#Isotropic Dependencies
from half_edge import HalfEdgeModel
from iso_rem import IsotropicRemesher


isMacOS = (platform.system() == "Darwin")

class AppWindow:

    MENU_OPEN = 1
    MENU_QUIT = 2

    def __init__(self, width, height):

        resource_path = gui.Application.instance.resource_path
        self.w_width = width
        self.w_height = height
        self.first_click = True

        #boilerplate - initialize window & scene
        self.window = gui.Application.instance.create_window("Test", width, height)
        self._scene = gui.SceneWidget()
        self._scene.scene = rendering.Open3DScene(self.window.renderer)

        #initializing menubar
        self._init_menubar()

        #basic layout
        self.window.set_on_layout(self._on_layout)
        self.window.add_child(self._scene)

        #set mouse and key callbacks
        self._scene.set_on_key(self._on_key_pressed)
        self._scene.set_on_mouse(self._on_mouse_pressed)
        
        #authorization for double presses
        self.event_auth = 0

        #geometry container for future reference
        self.geometry = None
        self.vertices = None
        self.vertex_normals = None
        self.triangles = None
        self.triangle_normals = None
        self.tree = None
        self.selected_vertex = None
        self.eigenvectors = None
        self.which = "mesh"
        self.current_eigenvector = 0

        #materials
        self.matlit = rendering.MaterialRecord()
        self.matlit.shader = "defaultLit"
        self.matunlit = rendering.MaterialRecord()
        self.matunlit.shader = "defaultUnlit"
        self.matline = rendering.MaterialRecord()
        self.matline.shader = "unlitLine"
        self.matline.line_width = 0.5

    def _init_menubar(self):

        if gui.Application.instance.menubar is None:

            if isMacOS:
                app_menu = gui.Menu()
                app_menu.add_item("Quit", AppWindow.MENU_QUIT)

            file_menu = gui.Menu()
            file_menu.add_item("Open...", AppWindow.MENU_OPEN)

            if not isMacOS:
                file_menu.add_separator()
                file_menu.add_item("Quit", AppWindow.MENU_QUIT)

            menu = gui.Menu()
            if isMacOS:
                # macOS will name the first menu item for the running application
                # (in our case, probably "Python"), regardless of what we call
                # it. This is the application menu, and it is where the
                # About..., Preferences..., and Quit menu items typically go.
                menu.add_menu("Example", app_menu)
                menu.add_menu("File", file_menu)
                # Don't include help menu unless it has something more than
                # About...
            else:
                menu.add_menu("File", file_menu)

            gui.Application.instance.menubar = menu

        # The menubar is global, but we need to connect the menu items to the
        # window, so that the window can call the appropriate function when the
        # menu item is activated.
        self.window.set_on_menu_item_activated(AppWindow.MENU_OPEN, self._on_menu_open)
        self.window.set_on_menu_item_activated(AppWindow.MENU_QUIT, self._on_menu_quit)

    def _on_menu_open(self):
        dlg = gui.FileDialog(gui.FileDialog.OPEN, "Choose file to load",
                             self.window.theme)
        dlg.add_filter(
            ".ply .stl .fbx .obj .off .gltf .glb",
            "Triangle mesh files (.ply, .stl, .fbx, .obj, .off, "
            ".gltf, .glb)")
        dlg.add_filter(
            ".xyz .xyzn .xyzrgb .ply .pcd .pts",
            "Point cloud files (.xyz, .xyzn, .xyzrgb, .ply, "
            ".pcd, .pts)")
        dlg.add_filter(".ply", "Polygon files (.ply)")
        dlg.add_filter(".stl", "Stereolithography files (.stl)")
        dlg.add_filter(".fbx", "Autodesk Filmbox files (.fbx)")
        dlg.add_filter(".obj", "Wavefront OBJ files (.obj)")
        dlg.add_filter(".off", "Object file format (.off)")
        dlg.add_filter(".gltf", "OpenGL transfer files (.gltf)")
        dlg.add_filter(".glb", "OpenGL binary transfer files (.glb)")
        dlg.add_filter(".xyz", "ASCII point cloud files (.xyz)")
        dlg.add_filter(".xyzn", "ASCII point cloud with normals (.xyzn)")
        dlg.add_filter(".xyzrgb",
                       "ASCII point cloud files with colors (.xyzrgb)")
        dlg.add_filter(".pcd", "Point Cloud Data files (.pcd)")
        dlg.add_filter(".pts", "3D Points files (.pts)")
        dlg.add_filter("", "All files")

        # A file dialog MUST define on_cancel and on_done functions
        dlg.set_on_cancel(self._on_file_dialog_cancel)
        dlg.set_on_done(self._on_load_dialog_done)
        self.window.show_dialog(dlg)
    
    def _on_file_dialog_cancel(self):
        self.window.close_dialog()
    
    def _on_load_dialog_done(self, filename):
        self.window.close_dialog()
        self.load(filename)

    def _on_menu_quit(self):
        gui.Application.instance.quit()

    def _preprocess(self, m):

        vertices, triangles = np.asarray(m.vertices), np.asarray(m.triangles)
        vertex_normals = np.asarray(m.compute_vertex_normals())
        triangle_normals = np.asarray(m.compute_triangle_normals())
        #centering
        vertices = vertices - vertices.mean(0)
        
        #unit_sphere_normalization
        norm = np.max((vertices * vertices).sum(-1))
        vertices = vertices / np.sqrt(norm)

        return o3d.geometry.TriangleMesh(o3d.utility.Vector3dVector(vertices), o3d.utility.Vector3iVector(triangles))

    def _find_match(self, query):

        if self.geometry is not None:

            ind = 0
            if self.tree is not None:
                _, ind, _ = self.tree.search_knn_vector_3d(query, 1)
                ind = int(np.asarray(ind)[0])
                self.selected_vertex = ind
                return self.vertices[ind]
                
            else:
                d = self.vertices - query
                d = np.argmin((d * d).sum(-1))
                self.selected_vertex = d
                return self.vertices[ind]

    def load(self, path = None):

        #clearing scene
        self._scene.scene.clear_geometry()
        if (path):
            #reading geometry type
            geometry_type = o3d.io.read_file_geometry_type(path)

            #checking the type of geometry
            if geometry_type & o3d.io.CONTAINS_TRIANGLES:
                self.geometry = o3d.io.read_triangle_model(path).meshes[0].mesh
                
        if self.geometry is None:
            print("[Info]", path, "appears to not be a triangle mesh")
            return 
        else:
            #preprocessing and setting geometry
            self.geometry = self._preprocess(self.geometry)
            self.line_geometry = o3d.geometry.LineSet.create_from_triangle_mesh(self.geometry)
            # self.line_geometry.paint_uniform_color([1, 0, 0])
            #setting vertex and triangle data for easy access
            self.vertices = np.asarray(self.geometry.vertices)
            self.triangles = np.asarray(self.geometry.triangles)
            self.vertex_normals = np.asarray(self.geometry.vertex_normals)
            self.triangle_normals = np.asarray(self.geometry.triangle_normals)
            #initializing kd-tree for quick searches
            self.tree = o3d.geometry.KDTreeFlann(self.geometry)

            #adding mesh to the scene and reconfiguring camera
            self._scene.scene.add_geometry("__model__", self.geometry, self.matlit)
            bounds = self._scene.scene.bounding_box
            self._scene.setup_camera(60, bounds, bounds.get_center())

    def _on_layout(self, layout_context):
        
        r = self.window.content_rect
        self._scene.frame = r

    def _on_key_pressed(self, event):

        print("key pressed: ", event.key)

        #number key -> ring neighborhood
        if (event.key <= 57 and event.key >= 48) and event.type == KeyEvent.Type.UP:
            self._show_k_ring(event.key-48)

        #C key - delta coordinates
        if event.key == 99 and event.type == KeyEvent.Type.UP:
            self._show_delta_coordinates()
            print("Delta Coordinates Visualization.")
            return gui.Widget.EventCallbackResult.HANDLED

        #S key - eigendecomposition and 
        elif event.key == 115 and event.type == KeyEvent.Type.UP:
            self._show_eigendecomposition()
            print("Eigenvisualization.")
            return gui.Widget.EventCallbackResult.HANDLED
        
        #V key - reset geometry and redraw scene
        elif event.key == 118 and event.type == KeyEvent.Type.UP:
            self._reset_gepometry()
            self._redraw_scene()
            return gui.Widget.EventCallbackResult.HANDLED
        
        #T key - toggle mesh or lineset
        elif event.key == 116 and event.type == KeyEvent.Type.UP:
            self.which = "line" if self.which == "mesh" else "mesh"
            self._reset_gepometry()
            self._redraw_scene()
            
            # self._scene.scene.add_geometry("__model_lined__", self.geometry, self.matline)
            print("mode = ", self.which)
            return gui.Widget.EventCallbackResult.HANDLED
        
        #R key - eigenvector visualization mode
        elif event.key == 114 and event.type == KeyEvent.Type.UP:
            self._calc_eigenvectors()
            print("eigenvectors calculated.")
            return gui.Widget.EventCallbackResult.HANDLED

        #left or bottom arrow keys - decrease eigenvector counter
        elif (event.key == 263 or event.key == 266) and event.type == KeyEvent.Type.UP:
            self.current_eigenvector = self.current_eigenvector -1 if self.current_eigenvector > 0 else 0
            print("current eigenvector: ", self.current_eigenvector)
            return gui.Widget.EventCallbackResult.HANDLED
        
        #right or up arrow keys - increase eigenvector counter
        elif (event.key == 264 or event.key == 265) and event.type == KeyEvent.Type.UP:
            self.current_eigenvector = self.current_eigenvector +1 if self.current_eigenvector < self.vertices.shape[0]-1 else self.vertices.shape[0]-1
            print("current eigenvector: ", self.current_eigenvector)
            return gui.Widget.EventCallbackResult.HANDLED

        #enter key - show eigenvector
        elif event.key == 10 and event.type == KeyEvent.Type.UP:
            self._show_eigenvector()
            return gui.Widget.EventCallbackResult.HANDLED
        
        ##MY ADDED BUTTONS
        # L for laplacian
        #K for taubin
        elif event.key == 108 and event.type == KeyEvent.Type.UP:
            for i in range(6):
                self._laplace_smoothing()
            print("laplace smoothing")
            return gui.Widget.EventCallbackResult.HANDLED

        elif event.key == 107 and event.type == KeyEvent.Type.UP:
            for i in range(6):
                self._taubin_smoothing()
            print("taubin smoothing")
            return gui.Widget.EventCallbackResult.HANDLED

        #B for Ball Pivot Reconstruction
        elif event.key == 98 and event.type == KeyEvent.Type.UP: 
            self._ball_pivot_reconstruction()
            print("Reconstructed using vertices and ball pivot method")
            return gui.Widget.EventCallbackResult.HANDLED
        
        #U for curvature mapping
        elif event.key == 117 and event.type == KeyEvent.Type.UP: 
            self._curvature_mapping()
            print("Colored mesh according to Curvature")
            return gui.Widget.EventCallbackResult.HANDLED
        
        #I for Isotropic Remeshing
        elif event.key == 105 and event.type == KeyEvent.Type.UP:
            self._isotropic_remeshing()
            
            print("Mesh just became meshier.")
            return gui.Widget.EventCallbackResult.HANDLED
        
        #G for Multiresolutional Reeb Graph Analysis
        elif event.key == 103 and event.type == KeyEvent.Type.UP:
            self._topological_analysis()
            print("Extracted topological data.")
            return gui.Widget.EventCallbackResult.HANDLED
        

        return gui.Widget.EventCallbackResult.IGNORED
    

    def _on_mouse_pressed(self, event):

        # interfere with manipulating the scene.
        if event.type == gui.MouseEvent.Type.BUTTON_DOWN and event.is_modifier_down(gui.KeyModifier.CTRL):

            def depth_callback(depth_image):
                # Coordinates are expressed in absolute coordinates of the
                # window, but to dereference the image correctly we need them
                # relative to the origin of the widget. Note that even if the
                # scene widget is the only thing in the window, if a menubar
                # exists it also takes up space in the window (except on macOS).
                x = event.x - self._scene.frame.x
                y = event.y - self._scene.frame.y
                # Note that np.asarray() reverses the axes.
                depth = np.asarray(depth_image)[y, x]

                if depth == 1.0:  # clicked on nothing (i.e. the far plane)
                    self.selected_vertex = None
                else:
                    #finding point on the mesh
                    world = self._scene.scene.camera.unproject(
                        x, y, depth, self._scene.frame.width,
                        self._scene.frame.height)
                    
                    #finding closest mesh vertex
                    match = self._find_match(world)

                    #adding a sphere
                    sphere = o3d.geometry.TriangleMesh.create_sphere(radius=0.01).paint_uniform_color(np.array([1,0,0])).translate(match)
                    self._scene.scene.remove_geometry("__point__")
                    self._scene.scene.add_geometry("__point__", sphere, self.matlit)


            self._scene.scene.scene.render_to_depth_image(depth_callback)
            return gui.Widget.EventCallbackResult.HANDLED
    
        elif event.type == MouseEvent.Type.BUTTON_DOWN:
            return gui.Widget.EventCallbackResult.HANDLED
        
        
        
        return gui.Widget.EventCallbackResult.HANDLED

    def _set_projection(self):
        self._scene.scene.camera.set_projection(
            Camera.Projection(1), -2.0, 2.0, -2.0, 2.0, 0.1, 100.0
        )

    def _reset_gepometry(self):

        self.geometry = o3d.geometry.TriangleMesh(
            o3d.utility.Vector3dVector(self.vertices),
            o3d.utility.Vector3iVector(self.triangles)
        )

        self.line_geometry = o3d.geometry.LineSet.create_from_triangle_mesh(self.geometry)
        self.line_geometry.paint_uniform_color([0, 0, 0])

    def _redraw_scene(self):

        #clearing scene
        self._scene.scene.clear_geometry()
        
        #if line mode then draw lineset
        if self.which == "line":
            self._scene.scene.add_geometry("__model__", self.line_geometry, self.matline)
        elif self.which == "mesh":
            self._scene.scene.add_geometry("__model__", self.geometry, self.matlit)

    def _show_delta_coordinates(self):

        #calculate delta coordinates
        delta = U.delta_coordinates(self.vertices, self.triangles, use_laplacian=True)

        #calculating norm of delta vector
        norm = np.sqrt((delta * delta).sum(-1))

        #linear transformation
        norm = (norm - norm.min()) / (norm.max() - norm.min())

        #coloring the mesh
        colors = U.sample_colormap(norm)
        # colors = np.zeros_like(self.vertices)
        # colors[:,0] = norm

        self.geometry.vertex_colors = o3d.utility.Vector3dVector(colors)

        self._redraw_scene()

    def _show_k_ring(self, k):
        self.selected_vertex = 1
        if self.selected_vertex is not None:
            print(f"finding {k}-ring neighbors")
            neighbors = U.k_ring_adjacency(self.selected_vertex, self.triangles, k)
            # neighbors = U.k_ring_recursive(self.selected_vertex, self.triangles, k)

            colors = np.zeros_like(self.vertices)
            colors[neighbors, 0] = 1

            self.geometry.vertex_colors = o3d.utility.Vector3dVector(colors)

            self._redraw_scene()
        else:
            print("No vertex selected.")

    def _show_eigendecomposition(self):

        if self.geometry is not None:
            #constants
            num_components = self.vertices.shape[0]
            keep_percent = 0.1
            keep_components = int(num_components * keep_percent)
            discard_components = int(num_components * (1-keep_percent))

            #calculating the graph laplacian
            L = U.graph_laplacian(self.triangles).astype(np.float64)

            #eigen decomposition of symmetric matrix -> SM means return the smallest eigenvalues
            vals, vecs = eigsh(L, k=keep_components, which='SM')

            #forming the eigenvector matrix with only the significant components
            U_k = vecs#[:, 0:keep_components]
            V_filtered = U_k @ (U_k.T @ self.vertices)

            #setting the vertices to be the filtered ones
            self.geometry.vertices = o3d.utility.Vector3dVector(V_filtered)

            #redrawing to see the difference
            self._redraw_scene()

    def _calc_eigenvectors(self):

        if self.geometry is not None:
            
            #calculating the graph laplacian
            L = U.graph_laplacian(self.triangles).astype(np.float32).toarray()
            
            #performing eigendecomposition
            vals, vecs = eigh(L)

            #sorting according to eigenvalue
            sort_idx = np.argsort(vals)
            self.eigenvectors = vecs[:, sort_idx]

    def _show_eigenvector(self):
        
        if self.eigenvectors is not None:
            
            # colors = np.zeros_like(self.vertices)
            
            scalars = self.eigenvectors[self.current_eigenvector]
            scalars = (scalars - scalars.min()) / (scalars.max() - scalars.min())

            # colors[:,0] = scalars
            colors = U.sample_colormap(scalars)

            self.geometry.vertex_colors = o3d.utility.Vector3dVector(colors)
            self._redraw_scene()
    
    def _laplace_smoothing(self):
        if self.geometry is not None:
            l = 0.6
            vertices_smooth = U.laplacian_smoothing(self.triangles, self.vertices, l)

            self.geometry.vertices = o3d.utility.Vector3dVector(vertices_smooth)

            self._redraw_scene()
        
    def _taubin_smoothing(self):
        if self.geometry is not None:
            l = 0.6
            m = 0.3
            vertices_smooth = U.taubin_smoothing(self.triangles, self.vertices, l, m)

            self.geometry.vertices = o3d.utility.Vector3dVector(vertices_smooth)

            self._redraw_scene()

    def _ball_pivot_reconstruction(self):
        radii = [ 0.005, 0.008, 0.01, 0.014, 0.016, 0.005, 0.008, 0.01, 0.014, 0.016, 0.018, 0.024, 0.04]
        # radii = [0.006 * i for i in range(1, 20)]
        if self.geometry is not None:
            point_cloud = o3d.geometry.PointCloud(self.geometry.vertices)
            point_cloud.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius = 0.06,  max_nn=60))

            self.geometry = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(
            point_cloud, o3d.utility.DoubleVector(radii))

            self._redraw_scene()

    def _curvature_mapping(self):
        if self.geometry is not None:
            principal_curvature, principal_dir1, principal_dir2 = crv.GetCurvaturesAndDerivatives(self)
            mean_curv = 0.5 * (principal_curvature[0, :] + principal_curvature[1, :])
            colors = U.sample_colormap(mean_curv)
            self.geometry.vertex_colors = o3d.utility.Vector3dVector(colors)
            self._redraw_scene()

    def _isotropic_remeshing(self):
        if self.geometry is not None:
            model = HalfEdgeModel(self.vertices, self.triangles)
            remesher = IsotropicRemesher(model)
            
            L_min, L_mean, L_max = U.compute_len_stats(self.geometry.vertices)
            print(f'L min = {L_min:.2f} L mean = {L_mean:.2f} L max = {L_max:.2f}')

            L = 1*L_mean
            ITER = 1

            remesher.isotropic_remeshing(
                L, 
                iter=ITER, 
                explicit=True, 
                foldover=20,
                sliver=True
            )
            model.clean()
            _triangle_mesh = o3d.geometry.TriangleMesh(model.vertices, model.triangles)
            self.geometry = _triangle_mesh
            self.load()
            # self._reset_gepometry()
            # self._preprocess(_triangle_mesh)
            # self._scene.scene.add_geometry("__remeshed_model__", _triangle_mesh, self.matlit)
            # self._scene.scene.remove_geometry("__model__")
            # self._redraw_scene()
            # o3d.visualization.draw_geometries([o3d.geometry.PointCloud(self.geometry.vertices), _triangle_mesh], mesh_show_back_face=True, mesh_show_wireframe=True) 

    def _topological_analysis(self):
        if self.geometry is not None:
            S_area = self.geometry.get_surface_area()
            print("Mesh Surface Area: ", S_area)
            ##Shortcut edges are ommited, estimation error increased.
            # U.add_shortcut_edges(self.triangles, self.vertices)
            U.geodesic_dijkstra(self.vertices, self.triangles, S_area)
            print("Wait, that worked?")
        

def main():

    gui.Application.instance.initialize()

    w = AppWindow(1024, 768)

    if len(sys.argv) > 1:
        path = sys.argv[1]
        if os.path.exists(path):
            w.load(path)
        else:
            w.window.show_message_box("Error",
                                      "Could not open file '" + path + "'")

    # Run the event loop. This will not return until the last window is closed.
    gui.Application.instance.run()

if __name__ == "__main__":

    main()