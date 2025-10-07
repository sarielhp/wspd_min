# wpsd_min


This is a fork of my older code (currently private?) in Julia for computing 
[well-separated pairs decomposition](https://en.wikipedia.org/wiki/Well-separated_pair_decomposition).

Installation
------------
The easiest thing is to copy the scripts/julial script and put
it somewhere visible on the path. Then run (from the command line):
src/deps.jl. This should install all the needed packages in julia.


Use
---
running example/wspd_1_dim.jl should generate a figure in the out
directory. See the source for core more information. Note that it is
using Gurobi for the IP solving. Good luck installing that!


# Computing wspd for 1 dimension

Run 

> examples/wspd_1_dim.jl 1dim

From the "root" directory of the project.

This should create a figure out/1_dim.pdf

# System

I assume you use linux. If you use windows, install wsl, and use
linux. If you use mac you are on your own.


