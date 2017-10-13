function y = xlogx(x)
y = zeros(size(x));
for k = 1:length(x)
    if (x(k) > eps)
        y(k) = x(k)*log(x(k));
    end
end
