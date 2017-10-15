function res = numquad(f,a,b,N,mode)
% Numerical quadrature on $[a,b]$ by polynomial quadrature formula
% f -> function to be integrated (handle), must support vector arguments
% a,b -> integration interval $[a,b]$ (endpoints included)
% N -> Maximal degree of polynomial
% mode (equidistant, CC = Clenshaw-Curtis, Gauss) selects quadrature rule
% Required function gaussQuad computing nodes and weights of Gauss-Legendre
% quadrature rule on [-1,1]
if (nargin < 5), mode = 'equidistant'; end
res = [];

if strcmp(mode,'Gauss')
  for deg=1:N
    [gx,w] = gaussquad(deg);
    % Gauss points for \Blue{$[a,b]$}
    x = 0.5*(b-a)*gx+0.5*(a+b);
    y = feval(f,x);
    res = [res; deg, 0.5*(b-a)*dot(w,y)];
  end
elseif strcmp(mode,'CC')
  for deg=1:N
    [gx,w] = fclencurt(deg,-1,1);
    % CC nodes for \Blue{$[a,b]$}
    x = 0.5*(b-a)*gx+0.5*(a+b);
    y = feval(f,x);
    res = [res; deg, 0.5*(b-a)*dot(w,y)];
  end
else
  p = (N+1:-1:1);
  w = (b.^p - a.^p)./p;
  for deg=1:N
    % equidistant quadrature nodes  
    x = (a:(b-a)/deg:b);
    % ``Quick and dirty'' implementation through polynomial interpolation
    y = feval(f,x);
    poly = polyfit(x,y,deg);
    res = [res; deg, dot(w(N+1-deg:N+1),poly)];
  end
end
