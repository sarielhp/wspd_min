# WPSD.jl

Code in Julia for computing 
[well-separated pairs decomposition](https://en.wikipedia.org/wiki/Well-separated_pair_decomposition), and
approximating the diameter of a point set.

Installation
------------
Uses some Julia code from the Frechet distance Julia
[implementation](https://github.com/sarielhp/FrechetDist.jl). Currently
this is a bit of a hack - you need to create a symlink from
FrechetDist/src/ directory to *src/cg/*.

For fun, I also implemented the diameter approximation algorithm from 
[here](https://sarielhp.org/p/00/diam.html) which is surprisingly easy
(and fast too?).

This is currently a bit of hack. I would release it as an official
package if there is interest.

Results for approximate diameter 
-------


Random points spread uniformly on a sphere in various dimensions. The
preset approximation ratio is 1.1 (but the result are much better, so...):

[2d sphere](results/2d_sphere.md)

[3d sphere](results/3d_sphere.md)

[4d sphere](results/4d_sphere.md)

[5d sphere](results/5d_sphere.md)

[6d sphere](results/6d_sphere.md)
# wspd_min
