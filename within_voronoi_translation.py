#!/usr/bin/env python

"""
within_voronoi_translation.py: Move antennas uniformly within their voronoi cell.

Noise is often added to the GPS coordinates of antennas to hinter's an attacker
ability to link outside information to the released database. This code takes as
input a list of antennas location and moves them uniformly within their voronoi
cell. The noise added is proportional to the density of antennas in the region
while preserving the overall structure of the mesh. To avoid issues, all
proposed points have to fall within the convex hull formed by the antennas.

Use:
> import within_voronoi_translation as wvt
> wvt.generate_new_positions([(0.367, 0.491), (0.415, 0.289), (0.495, 0.851),...])

Test:
$ python within_voronoi_translation.py

Algorithm:
A delaunay triangulation is used to compute the neighbors of each centroids.
Points are then draw at random in the square bounding the circle who diameter is
equal to the distance between the centroid and its farthest neighbors.
Points are rejected until they fall into the voronoi cell and inside the convex
hull.
"""


import scipy.spatial
import random
import numpy as np


def __compute_distance(a, pos, positions):
  x1 = positions[a][0]
  y1 = positions[a][1]
  x2, y2 = pos[0], pos[1]
  return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)


def __compute_distance_centroids(a, b, positions):
  return __compute_distance(a, positions[b], positions)


def __compute_neighbors(positions):
  delaunay = scipy.spatial.Delaunay(positions)
  slices, delaunay_neighbors = delaunay.vertex_neighbor_vertices
  neighbors = []
  for node, pos in enumerate(positions):
    neighbors.append(list(delaunay_neighbors[slices[node]:slices[node + 1]]))
  return neighbors


def __compute_max_vertices(positions):
  voronoi = scipy.spatial.Voronoi(positions)
  neighboring_vertices = []
  for point, region in enumerate(voronoi.point_region):
    neighboring_vertices.append(max([__compute_distance(point, voronoi.vertices[pos], positions) for pos in voronoi.regions[region]]))
  return neighboring_vertices


def __outside_convexhull(point, initial_positions):
  """
  Return True if the point falls outside of the convex hull.
  """
  if set(scipy.spatial.ConvexHull(initial_positions + [point]).vertices) - set(scipy.spatial.ConvexHull(initial_positions).vertices):
    return True
  else:
    return False


def __compute_border_points(positions):
  """
  Return the set of points for which at least one of their voronoi vertice falls
  outside of the convex hull.
  """
  voronoi = scipy.spatial.Voronoi(positions)
  vertices_outside = set([-1])
  for i, vertice in enumerate(voronoi.vertices):
    if __outside_convexhull(vertice, initial_positions):
      vertices_outside.add(i)

  points_outside = set()
  for region_id, region in enumerate(voronoi.regions):
    if any(point in vertices_outside for point in region):
      points_outside.add(list(voronoi.point_region).index(region_id))
  return points_outside


def __compute_radiuses(positions, neighbors):
  radiuses = []
  for node, pos in enumerate(positions):
    radiuses.append(max([__compute_distance_centroids(node,i,positions) for i in neighbors[node]]))
  return radiuses


def __draw_point(node, positions, neighbors, radiuses):
  condition = True
  while condition:
    trans_x, trans_y = [(random.random() - .5) * radiuses[node] for i in range(2)]
    proposed_point = (positions[node][0] - trans_x, positions[node][1] - trans_y)
    condition = __compute_distance(node, proposed_point, positions) > min([__compute_distance(i, proposed_point, positions) for i in neighbors[node]])
    condition = condition + __outside_convexhull(proposed_point, positions)
  return proposed_point


def generate_new_positions(positions):
  neighbors = __compute_neighbors(positions)
  radiuses = [max(i) for i in zip(__compute_max_vertices(initial_positions), __compute_radiuses(positions, neighbors))]
  return [__draw_point(i, positions, neighbors, radiuses) for i in range(len(positions))]


if __name__ == '__main__':
  import matplotlib.pyplot as plt
  initial_positions = [(random.random(), random.random()) for i in range(50)]
  new_positions = generate_new_positions(initial_positions)
  fig = plt.figure(figsize=(10,9))
  scipy.spatial.voronoi_plot_2d(scipy.spatial.Voronoi(initial_positions), plt.gca())
  for i, pos in enumerate(initial_positions):
    plt.text(pos[0], pos[1], str(i))
  for i, pos in enumerate(new_positions):
    initial_pos = initial_positions[i]
    plt.plot([initial_pos[0], pos[0]], [initial_pos[1], pos[1]], 'k-')
    plt.plot(pos[0], pos[1], marker='o', color='r', ls='')
  for point in __compute_border_points(initial_positions):
    initial_pos = initial_positions[point]
    plt.plot(initial_pos[0], initial_pos[1], marker='o', color='g', ls='')
  hull = scipy.spatial.ConvexHull(initial_positions)
  for simplex in hull.simplices:
    list_x, list_y = zip(*[initial_positions[simplex[0]], initial_positions[simplex[1]]])
    plt.plot(list_x, list_y, 'b-')
  plt.show()

