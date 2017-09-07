function N1fx = evaluateN1(coordinates,elements,varargin)
%EVALUATEN1   Pointwise evaluation of the Newton potential N1.
%   EVALUATEN1 allows the pointwise evaluation (N1*fh)(z) of the Newton
%   potential N1*fh of a P0-function fh at points z on the boundary.
%
%   Usage: N1FX = EVALUATEN1(COORDINATES,ELEMENTS,VERTICES,VOLUMES,...
%                          FH,X,P2E,[,ETA])
%          N1FX = EVALUATEN1(COORDINATES,ELEMENTS,LAMBDA_H,X,P2E,[,ETA])
%
%   Let {E1,...,EK} be a partition of the boundary Gamma. COORDINATES gives
%   the coordinates for all nodes of the boundary partition. ELEMENTS 
%   specifies the boundary elements. Let {T1,...,TN} be a partition of the 
%   domain Omega. VERTICES gives the coordinates for all nodes of the 
%   partition, VOLUMES specifies the volume elements. The (N x 1)-vector FH
%   contains the elementwise values of the piecewise constantdensity fh, 
%   i.e., fh|_{Tj} = FH(j). Let X be an (M x 2)-matrix, where each row 
%   X(k,:) corresponds to some point x, where (N1*fh)(x) should be 
%   evaluated. Then, the function EVALUATEN1 returns an (M x 1)-vector N1FX
%   with N1FX(j) = N1*fh(X(j,:)). P2E is a (M x 1)-vector where P2E(j) = k 
%   is the index of the boundary element with X(j,:) in Ek. 
%
%   The optional parameter ETA is used for admissibilty criterions, cf. the
%   documentation of BUILDN, BUILDV and EVALUATKADJ.
% 
%   Note: The evaluation of (N1*fh)(z) at points z on the boundary is done 
%   via the formula
%      N1*fh(z) = [(-0.5+K')V^(-1)(N0f)](z)

% (C) 2009-2013 HILBERT-Team '09, '10, '12, '13
% support + bug report:  hilbert@asc.tuwien.ac.at
%
% Version: 3.1

if (nargin == 6)
    eta      = varargin{4};
elseif (nargin == 8)
    eta      = varargin{6};
else
    eta      = [];
end

if (nargin == 5 || nargin == 6)
    lambda_h = varargin{1};
    x        = varargin{2};
    p2e      = varargin{3};
elseif (nargin == 7 || nargin == 8)
    vertices = varargin{1};
    volumes  = varargin{2};
    fh       = varargin{3};
    x        = varargin{4};
    p2e      = varargin{5};

    N = buildN(coordinates,elements,vertices,volumes,eta);
    b = N*fh;
    clear N;

    V = buildV(coordinates,elements,eta);
    lambda_h = V\b;
    clear V;
else
    str = sprintf(strcat('Wrong number of input arguments! Either use\n',...
    '\tN1FX = evaluateN1(COORDINATES,ELEMENTS,VERTICES,VOLUMES,FH,X',...
    ')    or\n',...
    '\tN1FX = evaluateN1(COORDINATES,ELEMENTS,VERTICES,VOLUMES,FH,X',...
    ',ETA)\n'));
    error(str);
end

%*** Normalvector
cb1 = coordinates(elements(:,1),:);
cb2 = coordinates(elements(:,2),:);
nVec = [cb2(:,2)-cb1(:,2), cb1(:,1)-cb2(:,1)];
nVec = nVec./repmat(sqrt(sum(nVec.^2,2)),1,2);
nVec = nVec(p2e,:);

%*** Evaluate [(-0.5+K')lambda_h](x)
x_h = evaluateKadj(coordinates,elements,lambda_h,x,nVec,eta);
N1fx = -0.5*lambda_h(p2e)+x_h;
