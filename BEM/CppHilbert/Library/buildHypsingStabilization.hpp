function S = buildHypsingStabilization(coordinates,elements)
%BUILDHYPSINGSTABILIZATION   Stabilization matrix for the hypersingular IE.
%   S = BUILDHYPSINGSTABILIZATION(COORDINATES,ELEMENTS) assembles the 
%   stabilization matrix for the hypersingular IE and S1-elements.
%               
%   The hypersingular integral equation reads
%
%      W*u = (1/2-K')*phi
%
%   where phi is the known Neumann data for the solution u of
%
%      -Laplace(u) = f   in Omega
%            du/dn = phi on Gamma = boundary(Omega).
%
%   The corresponding stiffness matrix W of the hypersingular integral 
%   operator is usually built for S1 boundary elements. It therefore has a 
%   non-trivial kernel. This function computes the stabilization matrix
%
%      Sjk = (int_Gamma zetaj ds) (int_Gamma zetak ds)
%
%   which may be used to obtain a discrete system
%
%      (W+S)*u = RHS
%
%   which now allows for a uniquely determined numerical solution of the 
%   hypersingular integral equation.
%
%   COORDINATES gives the coordinates for all nodes on the boundary Gamma. 
%   ELEMENTS describes the boundary elements.

% (C) 2009-2013 HILBERT-Team '09, '10, '12, '13
% support + bug report:  hilbert@asc.tuwien.ac.at
%
% Version: 3.1

nE = size(elements,1);

%*** compute local mesh-size
h = sqrt(sum((coordinates(elements(:,1),:)...
    -coordinates(elements(:,2),:)).^2,2));

%*** build vector with entries c(j) = int_Gamma hatfunction(j) ds
c = 0.5*accumarray(reshape(elements,2*nE,1),[h;h]);

%*** build stabilization matrix
S = c*c';

