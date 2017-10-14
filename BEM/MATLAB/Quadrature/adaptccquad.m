function [I,n] = adaptccquad(f,a,b,atol,rtol)
% Adaptive Clenshaw-Curtis quadrature of f on [a,b]

% Minimum no. of quadrature points
n = 3;
% Maximal number of adaptive steps
maxsteps = 5;

[rx,w] = fclencurt(n,-1,1);
% CC nodes for \Blue{$[a,b]$}
x = 0.5*(b-a)*rx+0.5*(a+b);
y = feval(f,x);
% Quadrature formula: first evaluation
I = 0.5*(b-a)*dot(w,y); IH = I;

% Adaptive loop
for s=1:maxsteps
    n = 2*(n-1)+1;
    [rx,w] = fclencurt(n,-1,1);
    % CC nodes for \Blue{$[a,b]$}
    x = 0.5*(b-a)*rx+0.5*(a+b);
    % Inefficient implementation
    y = feval(f,x);
    I = 0.5*(b-a)*dot(w,y);
    est = abs(I - IH);
    fprintf('n = %i, Ih = %d, est = %i\n',n,I,est);
    % Stopping criterion
    if ((est < rtol*I) || (est < atol)), return; end
    IH = I;
end
