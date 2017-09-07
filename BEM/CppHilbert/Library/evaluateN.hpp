function Nf_x = evaluateN(varargin)
%EVALUATEN   Pointwise evaluation of the Newton potential.
%   NF_X = EVALUATEN(VERTICES,VOLUMES,FH,X) allows the pointwise evaluation
%   (N*fh)(x) of the Newton potential N*fh  of a P0-function fh at points x.
%
%   Let {T1,...,TN} be a partition of the domain Omega. VERTICES gives the 
%   coordinates for all nodes of the partition, VOLUMES specifies the volume
%   elements. The (N x 1)-vector FH contains the elementwise values of the 
%   piecewise constant density fh, i.e., fh|_{Tj} = FH(j). Let X be an
%   (M x 2)-matrix, where each row X(k,:) corresponds to some point x, where
%   (N*fh)(x) should be evaluated. Then, the function EVALUATEN returns an 
%   (M x 1)-vector
%      NF_X with NF_X(j) = N*fh(X(j,:)).

% (C) 2009-2013 HILBERT-Team '09, '10, '12, '13
% support + bug report:  hilbert@asc.tuwien.ac.at
%
% Version: 3.1

error('Please run MAKE in HILBERT root directory!');
