# Building graph for short characteristics method

The program builds graphs of cell traversal in the specified directions.

--- for input, the program accepts a setup file indicating the source files with data 
(grid, directions, 
output data (attention files are generated automatically. the number of files is equal to the number of directions),
and grid descriptor files

--- for the construction of working files (normals to cells, numbers of boundary cells, connection with neighboring cells, cells inside the border (if available)) the vtk library is needed

--- graph construction is available without using the vtk library

--- omp and mpi technologies are functionally connected