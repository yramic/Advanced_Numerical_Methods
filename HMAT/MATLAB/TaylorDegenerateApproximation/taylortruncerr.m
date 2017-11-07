function [ERR,X,Y] = taylortruncerr(q,x,y)
% Computes the approximation error of Taylor approximation of 1/(1+(x-y)^2) using
% Taylor expansion around 0. 
% q = number of summands, x,y = sampling points for error

[X,Y] = meshgrid(x,y);
ERR = 1./(1+(X-Y).^2);
for k=0:q
    S = (X-Y).^(2*k);
%    S = zeros(size(Y));
% $$$     for l=0:2*k
% $$$         S = S + binomial(2*k,l).*(X.^l).*((-Y).^(2*k-l));
% $$$     end
    ERR = ERR - ((-1)^k)*S;
end

% $$$ for l=0:2*q
% $$$   S = zeros(size(Y));
% $$$   for k=ceil(l/2):q
% $$$       S = S + ((-1)^k)*binomial(2*k,l)*((-Y).^(2*k-l));
% $$$   end
% $$$   ERR = ERR - (X.^l) .* S;
% $$$ end

