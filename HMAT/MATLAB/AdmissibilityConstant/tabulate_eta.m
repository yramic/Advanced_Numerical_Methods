% Matlab script for tabulating eta for pairs of boxes from exp:tprect
ret = [];
for k=0:5
    ret = [ret, [k;eta([0.55+k*0.05,0.75+k*0.05],[0.25-k*0.05,0.45-k*0.05])]];
end
ret,

ret = [];
am = 0.5*(sqrt(2)-1);
ap = 0.5*(sqrt(2)+1);
for k=0.05+0.04*(0:5)
    ret = [ret, [k;eta([k*am+0.5,k*ap+0.5],[-k*am+0.5,-ap*k+0.5])]];
end
ret,