function r = binomial(n,k)
% binomial coefficient 
r = 1;
for l=1:k
    r = r* (n-l+1)/l;
end
