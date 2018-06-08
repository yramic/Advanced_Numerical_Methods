function abeltest(L,plt)
% MATLAB script for testing implicit Euler convolution quadrature for the
% transfer function F(s) = s^(-1/2), which corresponds to the Abel integral
% operator

if (nargin < 2), plt = false; end
if (nargin < 1), L = 6; end

% Argument function
g = @(t) sqrt(t);
%g = @(t) (1-exp(-t));
% g = @(t) (1-exp(t));
%g = @(t) t;
% g = @(t) ones(size(t));
% Exact result
y = @(t) (sqrt(pi)*t/2);
% y = @(t) (2*(sqrt(t) - dawson(sqrt(t)))/sqrt(pi));
% y = @(t) (2*sqrt(t/pi) - exp(t).*erf(sqrt(t)));
% y = @(t) (4/3*t.*sqrt(t)/sqrt(pi));
% y = @(t) (2*sqrt(t)/sqrt(pi));

res = [];
for N=10*2.^(0:L)
   tau = 1.0/N; 
   tg = tau*(0:N);  % temporal grid
   w = cqrinvsqr(tau,N+1); % convolution quadrature weights
   gv = g(tg); % grid function g
   % Direct discrete convolution
   yh = zeros(size(tg));
   for l=1:N+1
       for k=1:l
          yh(l) = yh(l) + w(k)*gv(l-k+1); 
       end
   end
   % Compute error
   maxwerr = max(abs(yh-y(tg)));
   if (plt == true) 
      figure('name',sprintf('N=%d',N));
      plot(tg,yh,'r.',tg,y(tg),'b.');
      xlabel('t'); ylabel('y(t)/y_h(t)');
   end
   res = [res; N , maxwerr];
end

res,
figure('name','IE-CQ convergence');
loglog(1./res(:,1),res(:,2),'r-+',(1./res(:,1)),(1./res(:,1)*res(1,1)*res(1,2)*2),'k-');
xlabel('timestep \tau'); ylabel('err(\tau) (maximum norm)');
% title('g(t) = exp(-t) on [0,1]');
title('g(t) = sqrt(t) on [0,1]');
legend('err(\tau)','O(\tau)','location','best');

print -depsc2 'iecqcvg_sqrt.eps';
