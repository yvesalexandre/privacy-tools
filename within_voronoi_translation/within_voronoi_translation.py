#!/usr/bin/env python

"""
within_voronoi_translation.py: Move antennas uniformly within their voronoi cell.

Noise is often added to the GPS coordinates of antennas to hinter's an attacker
ability to link outside information to the released database. This code takes as
input a list of antennas location and moves them uniformly within their voronoi
cell and either the convex hull formed by the antennas or the polygon. The noise
added is proportional to the density of antennas in the region while preserving
the overall structure of the mesh.

Use:
> import within_voronoi_translation as wvt
> wvt.generate_new_positions([(0.367, 0.491), (0.415, 0.289), (0.495, 0.851),...])

Test:
$ python within_voronoi_translation.py
or
$ python within_voronoi_translation.py senegal

Algorithm:
Points are then draw at random in the square bounding the circle whose diameter
is equal to the maximum of the distance between the centroid its voronoi vertices
or the half-min distance with its neighbors for border points.
Points are rejected until they fall in the voronoi cell and either inside the
convex hull or the polygon.

Author: Yves-Alexandre de Montjoye
https://github.com/yvesalexandre/privacy-tools
"""


import scipy.spatial
import random
import numpy as np


def __compute_distance(a, pos, positions):
  """
  Return the distance between an antenna a and a point pos (tuple)
  """
  x1 = positions[a][0]
  y1 = positions[a][1]
  x2, y2 = pos[0], pos[1]
  return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)


def __compute_distance_centroids(a, b, positions):
  """
  Return the distance between two antennas a and b
  """
  return __compute_distance(a, positions[b], positions)


def __compute_border_points(positions):
  """
  Return the set of points for which at least one of their voronoi vertice falls
  outside of the convex hull.
  """
  voronoi = scipy.spatial.Voronoi(positions)
  vertices_outside = set([-1])
  for i, vertice in enumerate(voronoi.vertices):
    if not __in_convexhull(vertice, initial_positions):
      vertices_outside.add(i)

  points_outside = set()
  for region_id, region in enumerate(voronoi.regions):
    if any(point in vertices_outside for point in region):
      points_outside.add(list(voronoi.point_region).index(region_id))
  return points_outside


def __compute_max_radius(positions, neighbors):
  """
  Return a list of the maximum distances between an antenna and its voronoi
  vertices.

  Note: the half-min distance to its neighbors for border points
  """
  voronoi = scipy.spatial.Voronoi(positions)
  border_points = __compute_border_points(positions)
  radiuses = []
  for point, region in enumerate(voronoi.point_region):
    if point not in border_points:
      radiuses.append(max([__compute_distance(point, voronoi.vertices[pos], positions) for pos in voronoi.regions[region]]))
    else:
      radiuses.append(min([__compute_distance_centroids(point,i,positions) for i in neighbors[point]]) / 2)
  return radiuses


def __compute_neighbors(positions):
  """
  Return a list of the neighbors of every antenna.
  """
  delaunay = scipy.spatial.Delaunay(positions)
  slices, delaunay_neighbors = delaunay.vertex_neighbor_vertices
  neighbors = []
  for node, pos in enumerate(positions):
    neighbors.append(list(delaunay_neighbors[slices[node]:slices[node + 1]]))
  return neighbors


def __in_convexhull(point, initial_positions):
  """
  Return True if the point is inside the convex hull.
  """
  if set(scipy.spatial.ConvexHull(initial_positions + [point]).vertices) - set(scipy.spatial.ConvexHull(initial_positions).vertices):
    return False
  else:
    return True


def __draw_point(node, positions, neighbors, radiuses, polygon):
  """
  Return the new position of the antenna.
  """
  condition = True
  while condition:
    trans_x, trans_y = [(random.random() - .5) * radiuses[node] for i in range(2)]
    proposed_point = (positions[node][0] - trans_x, positions[node][1] - trans_y)
    in_voronoi = __compute_distance(node, proposed_point, positions) < min([__compute_distance(i, proposed_point, positions) for i in neighbors[node]])
    if in_voronoi:
      if polygon:
        if __in_polygon(proposed_point, polygon):
          return proposed_point
      else:
        if __in_convexhull(proposed_point, positions):
          return proposed_point
  return proposed_point


def generate_new_positions(positions, polygon=None):
  """
  Return the new position for all the antennas.

  Shapefile:
  polygon expects a lonlat polygon. Shapefiles can loaded in python using
  shapefile and can be converted to lonlat format using pyproj.transform and
  the appropriate projection (http://www.prj2epsg.org/search).
  """
  neighbors = __compute_neighbors(positions)
  radiuses = __compute_max_radius(positions, neighbors)
  output = []
  for point_id in range(len(positions)):
    output.append(__draw_point(point_id, positions, neighbors, radiuses, polygon))
  return output


def __in_polygon(point,poly):
  """
  Return whether a point is in a polygon.

  Ray-casting Algorithm
  Adapted from http://geospatialpython.com/2011/08/point-in-polygon-2-on-line.html
  """
  x, y = point
  # check if point is a vertex
  if (x,y) in poly:
    return True
  # check if point is on a boundary
  for i in range(len(poly)):
    p1 = None
    p2 = None
    if i == 0:
      p1 = poly[0]
      p2 = poly[1]
    else:
      p1 = poly[i - 1]
      p2 = poly[i]
    if p1[1] == p2[1] and p1[1] == y and x > min(p1[0], p2[0]) and x < max(p1[0], p2[0]):
      return True
  n = len(poly)
  inside = False
  p1x,p1y = poly[0]
  for i in range(n + 1):
    p2x,p2y = poly[i % n]
    if y > min(p1y,p2y):
      if y <= max(p1y,p2y):
        if x <= max(p1x,p2x):
          if p1y != p2y:
            xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
          if p1x == p2x or x <= xints:
            inside = not inside
    p1x,p1y = p2x,p2y
  return inside

if __name__ == '__main__':
  import sys
  import matplotlib.pyplot as plt
  if (len(sys.argv) > 1) and (sys.argv[1] == 'senegal'):
    main_type = False
  else:
    main_type = True
  if main_type:
    initial_positions = [(random.random(), random.random()) for i in range(350)]
    new_positions = generate_new_positions(initial_positions)
  else:
    import pyproj
    import shapefile
    sf = shapefile.Reader('senegal_shapefile/senegal.shp')
    region = sf.shapes()[0]
    polygon = [pyproj.transform(pyproj.Proj(init='epsg:32628'),pyproj.Proj(proj='latlong'), pts[0], pts[1]) for pts in region.points]
    initial_positions = []
    while len(initial_positions) < 350:
      point = [random.uniform(-18, -11), random.uniform(12,17)]
      if __in_polygon(point,polygon):
        initial_positions.append(point)
    new_positions = generate_new_positions(initial_positions, polygon)

  fig = plt.figure(figsize=(10,9))
  scipy.spatial.voronoi_plot_2d(scipy.spatial.Voronoi(initial_positions), plt.gca())
  for i, pos in enumerate(initial_positions):
    plt.text(pos[0], pos[1], str(i))
  for point in __compute_border_points(initial_positions):
    initial_pos = initial_positions[point]
    plt.plot(initial_pos[0], initial_pos[1], marker='o', color='g', ls='')
  if main_type:
    hull = scipy.spatial.ConvexHull(initial_positions)
    for simplex in hull.simplices:
      x, y = zip(*[initial_positions[simplex[0]], initial_positions[simplex[1]]])
      plt.plot(x, y, 'b-')
  else:
    list_x, list_y = zip(*polygon)
    plt.plot(list_x, list_y, 'b-')
  for i, pos in enumerate(new_positions):
    initial_pos = initial_positions[i]
    plt.plot([initial_pos[0], pos[0]], [initial_pos[1], pos[1]], 'k-')
    plt.plot(pos[0], pos[1], marker='o', color='r', ls='')
  plt.show()
