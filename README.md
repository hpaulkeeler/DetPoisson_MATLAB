# DetPoisson_MATLAB

A collection of MATLAB m. files used for the results in the paper  by Blaszczyszyn and Keeler[1]. In particular, they create and fit determinantally-thinned Poisson point process to dependently-thinned point processes such as Matern hard-core point processes; for details see[1]. The fitting method is based on that of Kulesza and Taskar[2].

Run the file DemoDetPoisson.m for a demonstration of a simulating a determinantally-thinned Poisson point process.

To reproduce the results in the paper by Blaszczyszyn and Keeler[1], first run SubsetGenerate.m, then SubsetDetPoissonFit.m, and finally SubsetDetPoissonGenerate.m. The first two files will create .mat files locally for storing variables and data. 

To reproduce the exact same results in the paper[2], set the random seed to one (ie add the line:rng(1)) in the files SubsetGenerate.m and SubsetDetPoissonFit.m. 

See comments in invidual .m files for more information. 

Author: H.P. Keeler, Inria/ENS, Paris, and University of Melbourne,
Melbourne, 2018.

References:
[1] Blaszczyszyn and Keeler, Determinantal thinning of point processes
with network learning applications, 2018.
[2] Kulesza and Taskar, "Determinantal point processes for machine learning",Now Publisers, 2012
