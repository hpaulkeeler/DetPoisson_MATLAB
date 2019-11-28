# DetPoisson_MATLAB

A determinantally-thinned (Poisson) point process is essentially a discrete determinantal point process whose underlying state space is a single realization of a (Poisson) point process defined on some (bounded) continuous space. I believe this is a new type of point process, originally proposed in the paper by Blaszczyszyn and Keeler[1]: 

https://arxiv.org/abs/1810.08672

An obvious question is whether a determinantally-thinned Poisson point process is *also* a determinantal point process? The answer, we believe, is no, but it's from obvious. 

## Demonstration

Run the file DemoDetPoisson.m for a demonstration of simulating/sampling a determinantally-thinned Poisson point process. 

## Reproducing results from the paper
Most of the remaining MATLAB files were used to obtain the results in the paper by Blaszczyszyn and Keeler[1]. In particular, they create and fit determinantally-thinned Poisson point process to dependently-thinned Poisson point processes such as Matern hard-core point processes; for details see[1]. The fitting (or supervised learning) method is based on that developed by Kulesza and Taskar[2].

To reproduce the results in the paper by Blaszczyszyn and Keeler[1], first run SubsetGenerate.m, then SubsetDetPoissonFit.m, and finally SubsetDetPoissonGenerate.m. The first two files will create .mat files (stored locally), which contain variable values (that is, simulation/sampling and fitting results). 

To reproduce the exact same results in the paper[1], set the random seed to one in the files SubsetGenerate.m and SubsetDetPoissonGenerate.m. That is, in these two files remove the comment from the line that reads:

%seedRand=1; rng(seedRand) 

See comments in the individual MATLAB .m files for more information. 

Author: H.P. Keeler, Inria/ENS, Paris, and University of Melbourne,
Melbourne, 2018.

## Other code repositories
I also wrote some of this code (for example DemoDetPoisson,m and TestDetPoisson) in R and in Python, which have a very similar structure; see  

https://github.com/hpaulkeeler/DetPoisson_Julia

https://github.com/hpaulkeeler/DetPoisson_Python

https://github.com/hpaulkeeler/DetPoisson_R 


## References

[1] Blaszczyszyn and Keeler, Determinantal thinning of point processes
with network learning applications, 2018.

https://arxiv.org/abs/1810.08672

[2] Kulesza and Taskar, "Determinantal point processes for machine learning",Now Publisers, 2012
