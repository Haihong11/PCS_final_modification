## Composite-Body Algorithm (CBA) of Discrete Cosserat dynamics model 


##This MATLAB code is aimed at implementing the composite-body algorithm of the Discrete Cosserat approach for trunk robot dynamics,
 so that we can design the controller based on the Cosserat dynamics model.


1. piecewise_driver:  We can run it to simulate the dynamics. The trunk robot is divided into three unequal sections and spacing among discs, 
and the geometrical model can be set according the actual model. 

2. piecewise_CBA:  it can be consist of  several parts, and detailed notes have been added in the original code.

3. Visualization: We can obtain clear configuration and motion  of robot by running it. 

4. 'vector_tilde' and 'vector_hat' represent the isomorphism between the twist vector representation and
 the matrix representation of the Lie algebra.

5. 'piecewise_expmap', 'piecewise_invAdjoint' and 'piecewise_ADJ' are exponential representation of SE(3) and its Adjoint representation.

6. 'matrix_coAdjoint', 'matrix_coadj', 'matrix_Adjoint' and 'matrix_adj' are standard representation.

## Cite

If you use the code for scientific publication purpose, please cite:   ×××××××××××××

## Author

The author of the algorithm and MATLAB implementation is: 

Haihong Li (haihong.li@inria.fr)

Inria Lille - Nord Europe and CRIStAL - Centre de Recherche en 
Informatique Signal et Automatique de Lille, 59650 Villeneuved’Ascq, France