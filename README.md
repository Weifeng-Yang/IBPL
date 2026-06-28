# IBPL
This package contains code for the sparse non-negative matrix factorization with $\ell_0$-constraints ($\ell_0$-SNMF) problem and the sparse non-negative CP tensor decomposition with $\ell_0$-constraints ($\ell_0$-SNCP) problem in our paper[<sup>1</sup>](#refer-id). 

## Matlab code
"L0SNMF" is the sparse non-negative matrix factorization with $\ell_0$-constraints ($\ell_0$-SNMF) function, and "L0SNCP" is the sparse non-negative CP tensor decomposition with $\ell_0$-constraints ($\ell_0$-SNCP) function. 

## How to run our code
A toy example explains how to use the these function. For "L0SNMF", just run the function '[main_Run_me.m](L0SNMF/main_Run_me.m)'. 

For "L0SNCP", before running it, first add the toolbox 'tensortoolbox'[<sup>2</sup>](#refer-id) (www.tensortoolbox.org) to the running path of matlab, and then run the function '[main_Run_me.m](L0SNCP/main_Run_me.m)'. 

## Data
This code has built-in the data mentioned in our paper[<sup>1</sup>](#refer-id). 

## Reference
<div id="refer-id"></div>
[1] An Inertial Block Proximal Linearized  Method with Adaptive Momentum for Nonconvex and Nonsmooth Optimization. 

[2] Brett W. Bader and Tamara G. Kolda. 2006. Algorithm 862: MATLAB tensor classes for fast algorithm prototyping. ACM Trans. Math. Softw. 32, 4 (December 2006), 635–653. https://doi.org/10.1145/1186785.1186794
