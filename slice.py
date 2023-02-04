import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from shapely.geometry import GeometryCollection, LineString, Point, Polygon
from shapely import polygonize, polygonize_full, get_parts
from stl import mesh

from util import *

LAYER_HEIGHT = 2

# Load the STL file into memory
your_mesh = mesh.Mesh.from_file('test-objects/3DBenchy.stl')

# Create a plot
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
        
# Get the height of the mesh    
z_min, z_max = stl_get_height(your_mesh)

# Create an array of layer heights
layers = np.arange(z_min, z_max, LAYER_HEIGHT)

# Define the slicing plane
plane_normal = np.array([0, 0, 1])

# Optimize
buckets = sort_tris(your_mesh, layers)


#layers = [0]
# Slice
for layer, z_slice in enumerate(layers):
    print(f"layer {layer}, {z_slice}")
    
    # Move the slicing plane to this specific Z height
    plane_point = np.array([0, 0, z_slice])
    
    # Get the polygon representation of the cross-section
    # TODO: Remov
    #cross_section = polygonize_full(slice_stl(buckets[layer], plane_normal, plane_point))
    cross_section = slice_stl(buckets[layer], plane_normal, plane_point)
    #print(cross_section)
    # Plot the polygon
    for geom in cross_section:
        #print(geom)
        x, y = geom.coords.xy
        ax.plot(x, y, z_slice)
    #    print(geom)
    #    x, y = geom.exterior.coords.xy
    #    ax.plot(x, y, z_slice)
    
set_axes_equal(ax)
plt.show()