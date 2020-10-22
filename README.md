# triinv
### Inversion of displacement and stress data using triangular dislocation elements in Matlab

Based on algorithms published in:

Meade, B. J. (2007), Algorithms for the calculation of exact displacements, strains, and stresses for triangular dislocation elements in a uniform elastic half space, *Computers and Geosciences*, **33**, 1064–1075, [doi:10.1016/j.cageo.2006.12.003](http://dx.doi.org/10.1016/j.cageo.2006.12.003).

__triinv__ uses [__tridisl__](https://github.com/jploveless/tridisl) as a submodule. To clone __triinv__, run

    $ git clone --recursive https://github.com/jploveless/triinv.git
    
__triinvx.m__ is the main function and requires a minimum of three input arguments: 
  
    u = triinvx(p, s, beta);
    
where __p__ is a structure containing information about the triangular dislocation elements, __s__ is a structure containing the constraining data, and __beta__ defines the weighting of regularization applied in the slip estimation. Slip is returned to the vector __u__. Information about these and optional input arguments can be found in the [wiki](https://github.com/jploveless/triinv/wiki).
