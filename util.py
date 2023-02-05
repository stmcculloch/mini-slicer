import numpy as np
import shapely
from shapely.geometry import (GeometryCollection, LineString, MultiLineString,
                              Point, Polygon)


def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

def slice_stl(tris, normals, plane_normal, plane_point):
    """
    Returns a list of points that represent the polygon of the cross-section
    of the mesh along the given plane
    
    tris : stl.mesh.Mesh
        triangles in the format of Mesh.vectors
    normals : stl.mesh.normals
        normals in the format of Mesh.normals
    plane_normal : np.array
        Normal vector of the slicing plane
    plane_point : np.array
        Point on the slicing plane
    """
    polygon_edges = []
    for i, (v1, v2, v3) in enumerate(tris):
        
        # Edges are now vectors
        edges = np.array([v2 - v1, v3 - v2, v1 - v3])
        verts = np.array([v1, v2, v3])
       
        # Find the intersections between the triangles and the plane
        intersections = []
        for j, edge in enumerate(edges):
            # Skip if the edge and plane are parallel
            # TODO: May increase runtime
            if np.dot(edge, plane_normal) == 0:
                continue
            # Find the intersection of the edge and the plane
            t = np.dot(plane_normal, plane_point - verts[j]) / np.dot(plane_normal, edge)
            if t >= 0 and t <= 1:
                intersections.append(verts[j] + t * edge)
        # If the "edge" has the start and end point at the same location, then it is just a point, so we can discard it   
        # E.g. (1.5 0, 1.5 0) is a point, not an edge    
        if intersections and not shapely.equals_exact(Point(intersections[0]), Point(intersections[1]), tolerance=1e-3): 
            new_edge = intersections
            new_edge_vector = np.array(new_edge[1]) - np.array(new_edge[0])
            cross_product = np.cross(new_edge_vector, normals[i])
            # TODO: May need to do normals[i] dot product with slicing plane to project normal onto plane before doing cross product
            # 
            if cross_product[2] < 0:
                new_edge.reverse()
            
            # Remove the Z values
            for i, edge in enumerate(new_edge):
                new_edge[i] = edge[:2]
            #new_edge = shapely.set_precision(LineString(new_edge), 0.01)
            new_edge = (LineString(new_edge))
            polygon_edges.append(new_edge)
    return shapely.set_precision(MultiLineString(polygon_edges), 0.01)

def stl_get_height(mesh):
    vertices = np.array(mesh.v0)
    # Find the minimum and maximum Z values
    min_z = np.min(vertices[:,2])
    max_z = np.max(vertices[:,2])
    print(f"Min Z = {min_z}, max Z = {max_z}")
    
    return min_z, max_z

def sort_tris(mesh, layers):
    """
    Sorts the triangles in the mesh into buckets based on the layers they intersect
        
    params
    mesh: stl.mesh.Mesh
        stl mesh to slice
    layers: np.array
        array of layer heights to slice at
        
    returns
    buckets: list of lists
        list of lists of triangles that intersect each layer"""
    tri_buckets = [[] for _ in range(len(layers))]
    normal_buckets = [[] for _ in range(len(layers))]
    num_tris = len(mesh.vectors)    
    
    # TODO: Optimize this    
    # For each triangle, find the layer it intersects    
    for i in range(num_tris):
        v1, v2, v3 = mesh.vectors[i]
        # Find the minimum and maximum Z value for the triangle
        triangle_zs = np.array([v1[2], v2[2], v3[2]])
        min_z = np.min(triangle_zs)
        max_z = np.max(triangle_zs)
        
        # Find the layers that the triangle intersects        
        for j, layer in enumerate(layers):
            
            # The triangle intersects the layer if the triangle's minimum Z is 
            # less than the layer's Z and maximum Z is greater than the layer's Z
            if layer >= min_z and layer <= max_z:
                tri_buckets[j].append((v1, v2, v3))
                normal_buckets[j].append(mesh.normals[i])
    return tri_buckets, normal_buckets