import time

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from shapely import get_parts, polygonize, polygonize_full
from shapely.geometry import (GeometryCollection, LineString, MultiLineString,
                              Point, Polygon)
from stl import mesh

from util import *

start_time = time.time()
LAYER_HEIGHT = 0.999
LINE_WIDTH = 1
WALLS = 2
PLOT_OUTPUT = True
# Load the STL file into memory
your_mesh = mesh.Mesh.from_file('test-objects/3DBenchy.stl')

# Create a plot
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
#ax = fig.add_subplot()
        
# Get the height of the mesh    
z_min, z_max = stl_get_height(your_mesh)

# Create an array of layer heights
layers = np.arange(0, z_max, LAYER_HEIGHT)
    
# Define the slicing plane
plane_normal = np.array([0, 0, 1])

# Optimize
tris, normals = sort_tris(your_mesh, layers)
time_to_optimize = time.time() - start_time

# Slice
for layer, z_slice in enumerate(layers):
    print(f"layer {layer}, z = {z_slice}")
    
    # Move the slicing plane to this specific Z height
    plane_point = np.array([0, 0, z_slice])
    
    # Get the polygon representation of the cross-section
    polygons = slice_stl(tris[layer], normals[layer], plane_normal, plane_point)
    for poly in polygons.geoms:
        x, y = poly.exterior.xy
        ax.plot(x, y, z_slice, color='red', linewidth=0.25)

# Timing how long the program takes
time_to_slice = time.time() - start_time - time_to_optimize
print(f"{time_to_optimize:.2f} seconds to optimize\n{time_to_slice:.2f} seconds to slice\n{time_to_slice / len(layers):.2f} seconds per layer\n")

if PLOT_OUTPUT: 
    set_axes_equal(ax)
    plt.show()