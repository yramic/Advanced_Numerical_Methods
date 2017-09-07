function Kadjx = evaluateKadj(varargin)
%EVALUATEKADJ   Pointwise evaluation of adjoint double-layer potential.
%   KADJPHI_X = EVALUATEKADJ(COORDINATES,ELEMENTS,PHIH,X,N_X [,ETA]) allows
%   the pointwise evaluation K'*phih(x) of the adjoint double-layer 
%   potential K'*phih of a P0-function phih at points x.
%
%   Let {E1,...,EN} be a partition of a boundary piece Gamma. COORDINATES 
%   gives the coordinates for all nodes of the partition, ELEMENTS specifies
%   the boundary elements. The (N x 1)-vector PHIH contains the elementwise
%   values of the piecewise constant density phih, i.e., 
%   phih|_{Ej} = PHIH(j). Let X be an (M x 2)-matrix, where each row X(j,:)
%   corresponds to some evaluation point of (K')phi(x) at the boundary Gamma
%   and let NX be an (M x 2)-matrix, where each row contains the normal 
%   vector of the element which contains the corresponding evaluation point.
%   Then, the function EVALUATEKADJ returns an (M x 1)-vector KADJG_X with 
%   KADJG_X(j) = (K*)phih(X(j,:)).
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
