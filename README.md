# 3D-DSMC
Three Dimensional Direct Simulation Monte Carlo Solver

-----Debuging compile command-----
ifort -i4 -r8 -O0 -g -debug all -traceback -check bounds -fp-model source -debug inline-debug-info -zero -assume nounderscore dsmc_check.f90 dsmc_initialization.f90 dsmc_drifting.f90 dsmc_boundaries.f90 dsmc_indexing.f90 dsmc_collision.f90 dsmc_sampling.f90 dsmc_output.f90 dsmc3D.f90

-ffree-line-length-none -fdefault-real-8 -fopenmp -Wall -fbounds-check -finit-local-zero -fdump-parse-tree -fdump-core -fbacktrace -fdefault-double-8 -fbackslash -O0 -fcray-pointer -Wno-lto-type-mismatch



-----Execution run compile command-----

ifort -i4 -r8 -O3 -heap-arrays -traceback -fp-model source -debug inline-debug-info -warn interfaces -zero -assume nounderscore dsmc_check.f90 dsmc_initialization.f90 dsmc_drifting.f90 dsmc_boundaries.f90 dsmc_indexing.f90 dsmc_collision.f90 dsmc_sampling.f90 dsmc_output.f90 dsmc3D.f90
