function logquaderrs(quad)
% Numerical quadrature of logarithms

if (nargin < 1), quad = 'Gauss'; end

N = 20; % max no of quadrature points 

% Integrand
f = @(x,a) log(x+a);
F = @(x,a) (x+a)*(log(x+a)-1);

a = 1.4; exact = F(1.0,a)-F(-1.0,a);
res1 = numquad(@(x) f(x,a),-1.0,1.0,N,quad);
err1 = abs(res1(:,2)-exact);

a = 1.2; exact = F(1.0,a)-F(-1.0,a);
res2 = numquad(@(x) f(x,a),-1.0,1.0,N,quad);
err2 = abs(res2(:,2)-exact);

a = 1.1; exact = F(1.0,a)-F(-1.0,a);
res3 = numquad(@(x) f(x,a),-1.0,1.0,N,quad);
err3 = abs(res3(:,2)-exact);

a = 1.05; exact = F(1.0,a)-F(-1.0,a);
res4 = numquad(@(x) f(x,a),-1.0,1.0,N,quad);
err4 = abs(res4(:,2)-exact);

figure('name','logquad: Gauss');
semilogy(res1(:,1),err1,'b+-',...
         res2(:,1),err2,'g+-',...
         res3(:,1),err3,'m+-',...
         res4(:,1),err4,'r+-');
title(sprintf('%s quadrature of t->log(t+\\alpha) on [-1,1]',quad));
xlabel('{\bf Number of quadrature nodes}','fontsize',12);
ylabel('{\bf |quadrature error|}','fontsize',12);
legend('\alpha=1.4','\alpha=1.2','\alpha=1.1','\alpha=1.05','location','best');
grid on;

print -depsc2 'logquaderrs.eps';

