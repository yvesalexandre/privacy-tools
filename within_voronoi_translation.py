#!/usr/bin/env python

"""
within_voronoi_translation.py: Move antennas uniformly within their voronoi cell.
"""


import scipy.spatial
import random
import numpy as np


def compute_distance(a, pos, positions):
  x1 = positions[a][0]
  y1 = positions[a][1]
  x2, y2 = pos[0], pos[1]
  return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)


def compute_distance_centroids(a, b, positions):
  return compute_distance(a, positions[b], positions)


def compute_neighbors(positions):
  delaunay = scipy.spatial.Delaunay(positions)
  slices, delaunay_neighbors = delaunay.vertex_neighbor_vertices
  neighbors = []
  for node, pos in enumerate(positions):
    neighbors.append(list(delaunay_neighbors[slices[node]:slices[node + 1]]))
  return neighbors


def compute_radiuses(positions, neighbors):
  radiuses = []
  for node, pos in enumerate(positions):
    radiuses.append(max([compute_distance_centroids(node,i,positions) for i in neighbors[node]]))
  return radiuses


def draw_point(node, positions, neighbors, radiuses):
  condition = True
  while condition:
    trans_x, trans_y = [(random.random() - .5) * radiuses[node] for i in range(2)]
    proposed_point = (positions[node][0] - trans_x, positions[node][1] - trans_y)
    condition = compute_distance(node, proposed_point, positions) > min([compute_distance(i, proposed_point, positions) for i in neighbors[node]])
  return proposed_point


def generate_new_positions(positions):
  neighbors = compute_neighbors(positions)
  radiuses = compute_radiuses(positions, neighbors)
  return [draw_point(i, positions, neighbors, radiuses) for i in range(len(positions))]


if __name__ == '__main__':
  import matplotlib.pyplot as plt
  initial_positions = [(random.random(), random.random()) for i in range(50)]
  new_positions = generate_new_positions(initial_positions)
  fig = plt.figure(figsize=(10,9))
  scipy.spatial.voronoi_plot_2d(scipy.spatial.Voronoi(initial_positions), plt.gca())
  for i, pos in enumerate(initial_positions):
    plt.text(pos[0], pos[1], str(i))
  for i, pos in enumerate(new_positions):
    plt.text(pos[0], pos[1], str(i), color='r')
    plt.plot(pos[0], pos[1], marker='o', color='r', ls='')
  plt.show()


