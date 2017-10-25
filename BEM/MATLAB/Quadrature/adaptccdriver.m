function adaptccdriver()
% Numerical experiment with adaptice CC quadrature

% Tolerances
rtol = 1E-2;
atol = 1E-6;

% Shift parameter
l = 1;
aset = 1.05:0.05:2;
for a = aset
    f = @(x) log(x+a);
    F = @(x) (x+a)*(log(x+a)-1);
    exact = F(1) - F(-1);
    [I(l),n(l)] = adaptccquad(f,-1,1,atol,rtol);
    err(l) = abs(I(l)-exact);
    fprintf('n = %d, error = %d\n',n(l),err(l));
    l = l+1;
end

figure('name','adaptcc');
plot(aset,n,'b*');
xlabel('\alpha');
ylabel('no. of f-evaluations (*)');
title('Adaptive C.-C. quadrature of log(t+\alpha)');
yyaxis right;
semilogy(aset,err,'r+');
ylabel('quadrature error (+)');

print -depsc2 'adaptcc.eps';
