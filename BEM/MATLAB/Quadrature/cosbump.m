function y = cosbump(x)
y = zeros(size(x))
for k =1:length(x)
    z = 4*(x(k)-0.5);
    if (abs(z) < 1)
        y(k) = cos(z)^2;
    end
end
