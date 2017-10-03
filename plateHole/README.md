## Stress analysis of a plate with a hole

It's from the OpenFOAM User Guide version 5.0 Tutorial 2.2.

### Main features:

- **edges** keyword to define non-straight edges.

- blocks do not all have the same orientation, take care to math the directions for connecting blocks.

- **symmetryPlane** type patches.

- **grad(U)** gradient discretisation scheme using **leastSquares** method.

- solidDisplacementFoam solver (better study the source code of this solver).

### Issues

- Why stress boundary condition is wrote to *0/D* file? Then how to define displacement b.c.?

- Post-processing: postProcessing -func "components(sigma)" command and also
postProcessing -func "singleGraph", how to write such post-processing functions?
