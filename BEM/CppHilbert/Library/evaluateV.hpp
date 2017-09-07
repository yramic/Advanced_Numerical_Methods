function Vphi_x = evaluateV(varargin)
%EVALUATEV   Pointwise evaluation of the simple-layer potential.
%   VPHI_X = EVALUATEV(COORDINATES,ELEMENTS,PHIH,X [,ETA]) allows the 
%   pointwise evaluation V*phih(x) of the simple-layer potential V*phih of 
%   a P0-function phih at points x.
%
%   Let {E1,...,EN} be a partition of a boundary piece Gamma. COORDINATES 
%   gives the coordinates for all nodes of the partition, ELEMENTS specifies
%   the boundary elements. The (N x 1)-vector PHIH contains the elementwise
%   values of the piecewise constant density phih, i.e., 
%   phih|_{Ej} = PHIH(j). Let X be an (M x 2)-matrix, where each row X(j,:)
%   corresponds to some evaluation point of V*phih(x) within the domain 
%   Omega or its boundary Gamma. Then, the function EVALUATEV returns an
%   (M x 1)-vector VPHI_X with VPHI_X(j) = V*phih(X(j,:)).
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

