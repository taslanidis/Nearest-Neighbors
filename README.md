# Nearest Neighbors
Hashing and Searching for Vectors and Curves in C++<br>
Implemented different Algorithms for searching nearest neighbors on vectors and curves. Hashing, Grids, Projections.

## Vectors
- LSH
- HyperCube


### [LSH:](https://github.com/Fanarosss/Nearest_Neighbors/tree/master/src/LSH)
Hashing Vectors with amplified hash functions and use Hash Table for storage.<br>
Get approximate nearest neighbor - Or R-nearest neighbors, with M1 metric.

### [HyperCube:](https://github.com/Fanarosss/Nearest_Neighbors/tree/master/src/Hypercube)
Hashing Vectors with amplified hash functions and use Hyper Cube for storage.
Get approximate nearest neighbor Or R-nearest neighbors, with M1 metric


## Curves
- Grid LSH
- Projections

### [Grid LSH:](https://github.com/Fanarosss/Nearest_Neighbors/tree/master/src/Grid_Curves)
Curve Vectorization with Grid, then calling LSH for vectors on the vectorized curves. (Poor method)

### [Projections:](https://github.com/Fanarosss/Nearest_Neighbors/blob/master/src/Projections_Curves/RandomProjections.cpp)
Find all relevant traversals, and project the nearest traversals and use LSH on vectorized traversals, instead of vectorizing Curves. (Better method for curves)

### General Use functions for I\O and searching:
[Helper Functions](https://github.com/Fanarosss/Nearest_Neighbors/blob/master/src/Helper_Functions.cpp)

### Usage:
- Compilation<br>
*make all*
- Instructions for usage <br>
*./build/lsh.x*<br>
*./build/hypercube.x*

### Collaborators
[Konstantinos Athinaios](https://github.com/KostasA97)
<br>
[Theofanis Aslanidis](https://github.com/Fanarosss)

