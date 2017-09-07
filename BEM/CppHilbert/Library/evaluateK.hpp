function Kg_z = evaluateK(varargin)
%EVALUATEK   Pointwise evaluation of the double-layer potential.
%   KG_X = EVALUATEK(COORDINATES,ELEMENTS,GH,X [,ETA]) allows the pointwise
%   evaluation K*gh(z) of the double-layer potential K*gh of an S1-function
%   gh at a list of evaluation points x.
%
%   Let {E1,...,EN} be a partition of a boundary piece Gamma with nodes 
%   z1,...,zn. Note that there holds n=N for a closed boundary only, whereas
%   N<n in general. COORDINATES gives the coordinates for all nodes of the 
%   partition, ELEMENTS specifies the boundary elements. The (n x 1)-vector
%   GH contains the nodal values of the S1-density gh, i.e., gh(zj) = GH(j).
%   Let X be an (M x 2)-matrix, where each row X(j,:) corresponds to some 
%   evaluation point of K*gh(z). Then, the function EVALUATEK returns an 
%   (M x 1)-vector
%      KG_X with KG_X(j) = K*gh(X(j,:)). 
%
%   ETA is an optional parameter that defines whether an evaluation point 
%   XJ = X(J,:) and an element EK = ELEMENTS(K,:) are admissable, i.e. 
%   whether there holds
%      diam(EK) < ETA*dist(XJ, EK),
%   where diam(EK) denotes the diameter of the element EK and dist(XJ,EK) 
%   denotes the distance between XJ and EK.
%
%   For each evaluation point, the integral is written as a sum of 
%   element-wise integrals. If an evaluation point and an element are not 
%   admissible, we compute the integral by using analytic formulas and 
%   otherwise, we use Gaussian quadrature. This leads to increased numerical
%   stability.
%
%   In case that ETA is omitted, a sane default value is choosen. In case 
%   that ETA=0, all integrals are evaluated analytically.

% (C) 2009-2013 HILBERT-Team '09, '10, '12, '13
% support + bug report:  hilbert@asc.tuwien.ac.at
%
% Version: 3.1

error('Please run MAKE in HILBERT root directory!');
