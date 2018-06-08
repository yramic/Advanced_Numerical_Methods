function w = cqrinvsqr(tau,n)
% Returns the first n convolution quadrature weights for the transfer
% function F(s) = s^(-1/2) for timestep tau.
w = ones(1,n);
for l=2:n
   w(l) = w(l-1)*(l-1.5)/(l-1);
end
w = sqrt(tau)*w;
