function y = tentfn(x)
y = zeros(size(x))
for k =1:length(x)
    z = 4*(x(k)-0.5);
    if (abs(z) < 1)
        y(k) = 1-abs(z);
    end
end
