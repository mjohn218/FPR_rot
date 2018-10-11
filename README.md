# FPR_rot
Rotational rigid-body reaction-diffusion code

Performs reaction-diffusion dynamics for multi-site rigid molecules in 3D and 2D. This specific code places bound sites at contact given their current orientation. It is not yet generalized to take in user-defined orientations for 'snapping' molecules into place after association. It will do this snapping for the clathrin assembly problem simulated in the paper.  


to compile:
enter src directory
make

requirements:
g++  (or some other c++ compiler)
gsl (GNU Scientific library) must be installed.

