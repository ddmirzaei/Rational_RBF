%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Rational RBF-PU interpolation with polyharmonic splines %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-----------  By: Davoud Mirzaei, University of Isfahan, Iran, November 2020    -------------
--------------------------------------------------------------------------------------------
This file contains the Matlab code for the rational RBF-PU method of
 
"E. Farazandeh, D. Mirzaei, A Rational RBF Interpolation with Conditionally Positive 
 Definite Kernels, Adv. Comput. Math. (2021)"

****Just run 'main.m' to see the results

-Initially, this code provides the results for a 2D example but all functions except PolyMat 
   and ScatPoints2D work in all dimensions. 
-To generalize PolyMat for other dimensions you just need to switch over different cases for
   the MultiIndex vector. Or, you may use an integer partitioning algorithm to produce this 
   vector in arbitrary dimensions for a given polynomial order. 
-Of course, the similar functions ScatPoints1D and ScatPoints3D can be simply developed by
   the user. 
-To force the function Rational_RBF working for non-scalable kernels just call it with scaling
   value h=1. In this case, the new kernel should be added to the KerMat function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%