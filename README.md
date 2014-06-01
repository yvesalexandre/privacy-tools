privacy-tools, A set of privacy tools for metadata
==================================================

  * [within_voronoi_transation](https://github.com/yvesalexandre/privacy-tools/tree/master/within_voronoi_translation): Noise is often added to the GPS coordinates of antennas to hinter's an attacker ability to link outside information to the released database. This code takes as input a list of antennas location and moves them uniformly within their voronoi cell and the convex hull formed by the antennas. The noise added is proportional to the density of antennas in the region while preserving the overall structure of the mesh. [Example output](https://github.com/yvesalexandre/privacy-tools/blob/master/within_voronoi_translation/example.png) where blue and green points are the original antennas with their id, red points the new coordinates, and in blue the convex hull.

This set of tools by [Yves-Alexandre de Montjoye](http://deMontjoye.com) can be cited as 10.5281/zenodo.9900 and is under the MIT License.
