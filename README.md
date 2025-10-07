# wpsd_min

This is a fork of my older code (currently private?) in Julia for computing 
[well-separated pairs decomposition](https://en.wikipedia.org/wiki/Well-separated_pair_decomposition).

Installation
------------
Uses some Julia code from the Frechet distance Julia
[implementation](https://github.com/sarielhp/FrechetDist.jl). The
easiest thing is to copy the scripts/julial script and put somewhere
visible on the path. Then do (from the command line): src/deps.jl. 


Use
---
running example/wspd_1_dim should generate a figure in the out directory. See the source for core more information. Note that it is using Gurobi for the IP solving. Good luck installing that!


# wspd_min
