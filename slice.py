from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Point, Polygon, LineString, GeometryCollection
from stl import mesh

LAYER_HEIGHT = 1

# Load the STL file into memory
your_mesh = mesh.Mesh.from_file('test-objects/3DBenchy.stl')

# Create a plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

plt.tight_layout()

def slice_stl(mesh, plane_normal, plane_point, tol=0.0001):
    """
    Returns a list of points that represent the polygon of the cross-section
    of the mesh along the given plane
    
    mesh : stl.mesh.Mesh
        STL mesh to slice
    plane_normal : np.array
        Normal vector of the slicing plane
    plane_point : np.array
        Point on the slicing plane
    tol : float
        Tolerance for considering a point as on the plane
    """
    polygon_points = []
    culled = 0
    for i, (v1, v2, v3) in enumerate(mesh):
        # Check if all three vertices are on one side of the plane
        dists = np.dot(np.array([v1, v2, v3]) - plane_point, plane_normal)
        if np.all(dists > tol) or np.all(dists < -tol):
            culled += 1
            continue
        
        # If not, we need to find the intersections
        edges = np.array([v2 - v1, v3 - v2, v1 - v3])
        verts = np.array([v1, v2, v3])
        # Edges are now vectors
       
        intersections = []
        for j, edge in enumerate(edges):
            # Skip if the edge and plane are parallel
            if np.dot(edge, plane_normal) == 0:
                continue
            # Find the intersection of the edge and the plane
            t = np.dot(plane_normal, plane_point - verts[j]) / np.dot(plane_normal, edge)
            if t >= 0 and t <= 1:
                intersections.append(verts[j] + t * edge)
        polygon_points.extend(intersections)
    return polygon_points

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
    buckets = [[] for _ in range(len(layers))]
    num_tris = len(mesh.vectors)    
    
    # TODO: Optimize this    
    # For each triangle, find the layer it intersects    
    for i in range(num_tris):
        v1, v2, v3 = your_mesh.vectors[i]
    
        # Find the minimum and maximum Z value for the triangle
        triangle_zs = np.array([v1[2], v2[2], v3[2]])
        min_z = np.min(triangle_zs)
        max_z = np.max(triangle_zs)
        
        # Find the layers that the triangle intersects        
        for j, layer in enumerate(layers):
            
            # The triangle intersects the layer if the triangle's minimum Z is 
            # less than the layer's Z and maximum Z is greater than the layer's Z
            if layer >= min_z and layer <= max_z:
                buckets[j].append((v1, v2, v3))
    
    return buckets
        
# Get the height of the mesh    
z_min, z_max = stl_get_height(your_mesh)

# Create an array of layer heights
layers = np.arange(z_min, z_max, LAYER_HEIGHT)

# Define the slicing plane
plane_normal = np.array([0, 0, 1])

# Optimize
buckets = sort_tris(your_mesh, layers)

# Slice
for layer, z_slice in enumerate(layers):
    print(f"layer {layer}, {z_slice}")
    
    # Define the slicing plane
    plane_point = np.array([0, 0, z_slice])
    
    # Get the polygon representation of the cross-section
    polygon = slice_stl(buckets[layer], plane_normal, plane_point)
    
    # Skip if there are no points
    if(len(polygon) == 0):
        continue
    
    # Plot the polygon
    x, y, z = zip(*polygon)
    ax.plot(x, y, z_slice)
        
plt.show()