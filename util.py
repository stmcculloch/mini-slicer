import matplotlib.pyplot as plt
import numpy as np
import shapely
from shapely.geometry import (GeometryCollection, LineString, MultiLineString,
                              MultiPolygon, Point, Polygon)


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
        
        
        edges = np.array([v2 - v1, v3 - v2, v1 - v3]) # 3 vectors
        verts = np.array([v1, v2, v3]) # 3 points
        new_edge = []
        
        ## Need to deal with all cases of intersection:
        ## 1. All 3 vertices on plane
        ## 2. 2/3 vertices on plane
        ## 3. middle vertex on plane
        ## 4. top/bottom vertex on plane
        
        for j, edge in enumerate(edges):
            
            # Skip edge if the edge and plane are parallel
            if np.dot(plane_normal, edge) == 0:
                continue
            
            # Find the intersection of the edge and the plane
            t = np.dot(plane_normal, plane_point - verts[j]) / np.dot(plane_normal, edge)
            if 0 <= t <= 1:
                new_edge.append(verts[j] + t * edge)
        if(len(new_edge)) > 2:
                print("fuck")
                print(new_edge)        
        # If both points of intersections are too close together then ignore it
        if new_edge and not shapely.equals_exact(Point(new_edge[0]), Point(new_edge[1]), tolerance=1e-4): 
            
            new_edge_vector = np.array(new_edge[1]) - np.array(new_edge[0])
            cross_product = np.cross(new_edge_vector, normals[i])
            # TODO: May need to do normals[i] dot product with slicing plane to project normal onto plane before doing cross product
            # 
            if cross_product[2] < 0:
                new_edge.reverse()
            # Remove the Z values
            new_edge = [edge[:2] for edge in new_edge]
            
                
            # Add the edge to the list of edges    
            polygon_edges.append(LineString(new_edge))
    
    # At this point, we have a list of LineStrings that represent the edges of the cross-section:
    # [<LINESTRING (-1   -1   , -1    0.8)>, 
    #  <LINESTRING (-1    0.8 , -1    1  )>,
    #  <LINESTRING (-1    1   ,  0.8  1  )>, 
    #  <LINESTRING ( 0.8  1   ,  1    1  )>, 
    #  <LINESTRING ( 1    1   ,  1   -0.8)>, 
    #  <LINESTRING ( 1   -0.8 ,  1   -1  )>, 
    #  <LINESTRING (-0.8 -1   , -1   -1  )>, 
    #  <LINESTRING ( 1   -1   , -0.8 -1  )>]
    
    # Now we need to combine the edges into a single polygon
    # <LINESTRING (-1   -1   , -1    0.8, -1    1  ,  0.8  1  ,  1    1  ,  1   -0.8,  1   -1  , -0.8 -1  , -1   -1  )>
    # <POLYGON ((-1   -1   , -1    0.8, -1    1  ,  0.8  1  ,  1    1  ,  1   -0.8,  1   -1  , -0.8 -1  , -1   -1  ))>
    
    output_polys = []
    #for edge in polygon_edges:
    #    if len(list(edge.coords)) > 2:
    #        print("fuck\n\n\n\n\n\n\n\n")
    # Keep creating polygons until no more edges remain
    num_edges = len(polygon_edges)
    while(polygon_edges):
        # Create a new LINESTRING starting at the first edge: 
        # <LINESTRING (-1   -1   , -1    0.8)>
        
        
        poly = polygon_edges[0]
        
        # Then delete the first edge from the list
        polygon_edges.pop(0)
        poly_creation_failed = False
        # Repeat until the end of the LINESTRING is the same as the start of the LINESTRING
        i=0
        while(not shapely.equals_exact(Point(poly.coords[0]), Point(poly.coords[-1]), tolerance=1e-4)):
            i+=1
            
            if i % 10 == 0:
                print("Layer progress: ", round(100*(1-len(polygon_edges)/num_edges), -1), "%" , end='\r')
            # Look for the next edge that starts at the end of the current edge and add the second coordinate to the LINESTRING
            # <LINESTRING (-1   -1   , -1    0.8, -1    1  )>
            edge_found = False
            for j, next_edge in enumerate(polygon_edges):
                #print(len(list(next_edge.coords)))
                if shapely.equals_exact(Point(poly.coords[-1]), Point(next_edge.coords[0]), tolerance=1e-4):
                    poly = LineString(list(poly.coords) + list(next_edge.coords[1:]))
                    
                    # Then delete the next edge from the list
                    polygon_edges.pop(j)
                    edge_found = True
                    break
                if shapely.equals_exact(Point(poly.coords[-1]), Point(next_edge.coords[1]), tolerance=1e-4):
                    print(list(poly.coords[-2:]))
                    print(list(next_edge.coords[:]))
                    poly = LineString(list(poly.coords) + list(next_edge.coords[0:-1]))
                    
                    # Then delete the next edge from the list
                    polygon_edges.pop(j)
                    edge_found = True
                    break    
                
            if not edge_found:
                print("ERROR: Was not able to close polygon")
                print(poly)
                fig = plt.figure()
                ax = fig.add_subplot()
                ax.plot(*poly.xy, color='blue')
                for line in polygon_edges:
                    ax.plot(*line.xy, color='red')
                ax.axis('equal')
                plt.show()
                poly_creation_failed = True
                break
                
        #print(poly)        
        if poly_creation_failed:
            continue
        else:
            try:
                output_polys.append(Polygon(poly))
            except:
                print("ERROR: Was not able to create polygon")
                print(poly)
                fig = plt.figure()
                ax = fig.add_subplot()
                ax.plot(*poly.xy, color='blue')
                for line in polygon_edges:
                    ax.plot(*line.xy, color='red')
                ax.axis('equal')
                plt.show()
                continue
    print("                                      ", end='\r')    
        
    # Return the polygon edges as a MultiLineString
    # shapely.set_precision is used to remove floating point errors
    return MultiPolygon(output_polys)

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
        triangle_min_z = np.min(triangle_zs)
        triangle_max_z = np.max(triangle_zs)
        
        # Find the layers that the triangle intersects        
        for j, layer in enumerate(layers):
            
            # The triangle must intersect the plane if there is a point above and below the plane
            if triangle_min_z <= layer <= triangle_max_z:
                tri_buckets[j].append((v1, v2, v3))
                normal_buckets[j].append(mesh.normals[i])
                
    return tri_buckets, normal_buckets