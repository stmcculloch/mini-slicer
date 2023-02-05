import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from shapely import get_parts, polygonize, polygonize_full
from shapely.geometry import (GeometryCollection, LineString, MultiLineString,
                              Point, Polygon)
from stl import mesh

from util import *

LAYER_HEIGHT = 0.4
WALLS = 2
PLOT_OUTPUT = True
# Load the STL file into memory
your_mesh = mesh.Mesh.from_file('test-objects/Organizer.stl')

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

# Slice
for layer, z_slice in enumerate(layers):
    print(f"layer {layer}, {z_slice}")
    
    # Move the slicing plane to this specific Z height
    plane_point = np.array([0, 0, z_slice])
    
    # Get the polygon representation of the cross-section
    #cross_section = polygonize_full(slice_stl(buckets[layer], plane_normal, plane_point))
    polygons, cuts, dangles, invalids = polygonize_full(slice_stl(tris[layer], normals[layer], plane_normal, plane_point).geoms)

    wall_polys = []
    for i in range(WALLS):
        # TODO: positive buffer for holes
        wall_polys.append(polygons.buffer(-i/4))
    
    wall_polys = GeometryCollection(wall_polys)
    
    # Plot the polygon
    if PLOT_OUTPUT:   
        if wall_polys.geom_type == 'Polygon':
            x, y = polygons.coords.xy
            ax.plot(x, y, z_slice, color='green')
        for dangle in dangles.geoms:
            x, y = dangle.coords.xy
            ax.plot(x, y, z_slice, color='gray')
        else:    
            outer_wall = True
            for geom in wall_polys.geoms:
                if geom.geom_type == 'Polygon':
                    x, y = geom.exterior.coords.xy
                    if outer_wall:
                        ax.plot(x, y, z_slice, color='red', linewidth=0.25)
                        outer_wall = False
                    else:
                        ax.plot(x, y, z_slice, color='green', linewidth=0.25)
                else:    
                    for geom1 in geom.geoms:
                        x, y = geom1.exterior.coords.xy
                        if outer_wall:
                            ax.plot(x, y, z_slice, color='red', linewidth=0.25)
                            outer_wall = False
                        else:
                            ax.plot(x, y, z_slice, color='green', linewidth=0.25)
if PLOT_OUTPUT: 
    set_axes_equal(ax)
    plt.show()