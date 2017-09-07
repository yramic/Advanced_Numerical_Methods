function Wx = evaluateW(varargin)
%EVALUATEW   Pointwise evaluation of the hypersingular integral operator.
%   WG_X = EVALUATEW(COORDINATES,ELEMENTS,GH,X,N_X [,ETA]) allows the 
%   pointwise evaluation of W*gh(x) of the hypersingular integral operator W
%   for an S1-function gh at a list of points x.
%
%   Let {E1,...,EN} be a partition of a boundary piece Gamma. COORDINATES 
%   gives the coordinates for all nodes of the partition, ELEMENTS specifies
%   the boundary elements. The (N x 1)-vector GH contains the nodal values 
%   of the S1-density gh, i.e., gh(zj) = GH(j). Let X be an (M x 2)-matrix,
%   where each row X(j,:) corresponds to some evaluation point of W*gh(z) 
%   and let NX be an (M x 2)-matrix, where each row contains the normal 
%   vector of the element which contains the corresponding evaluation point.
%   Please note that each evaluation point has to be contained in exactly 
%   one boundary element, i.e. X(j,:) must not be equal to some node 
%   COORDINATES(i,:) for any i. Then, the function EVALUATEW returns an 
%   (M x 1)-vector
%      WG_X with WG_X(j) = W*gh(X(j,:)).
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

