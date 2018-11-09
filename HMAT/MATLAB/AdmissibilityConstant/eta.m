function e = eta(boxx,boxy)
% Computes the admissibility constant for two 1D intervals specified
% by 2-vectors.
if (length(boxx) ~= 2), error('boxx must be 2-vector'); end
if (length(boxy) ~= 2), error('boxy must be 2-vector'); end

if (boxx(1) > boxy(2)), tmp = boxx; boxx = boxy; boxy = tmp; end
if (boxx(2) > boxy(1)), error('Overlapping boxes'); end

dist = boxy(1)-boxx(2);
dx = boxx(2)-boxx(1);
dy = boxy(2)-boxy(1);

e = max(dx,dy)/dist;
    