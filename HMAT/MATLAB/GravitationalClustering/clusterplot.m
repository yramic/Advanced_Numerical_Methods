function clusterplot(bbox,X,masses)
% plot bounding box passed in the 4-vector bbox and the center f gravity of
% the stars contained in it. X is a n x 2 matrix containing the star
% positions. masses is a row vector of star masses
% x_min = bbox(1), x_max = bbox(2), y_min = bbox(3), y_max = bbox(4)
plot([bbox(1),bbox(2),bbox(2),bbox(1),bbox(1)],...
     [bbox(3),bbox(3),bbox(4),bbox(4),bbox(3)],'m-');
 idx = find((X(1,:)>= bbox(1)) .* (X(1,:)<=bbox(2)) .* ...
            (X(2,:) >= bbox(3)) .* (X(2,:)<=bbox(4)));
 tm = sum(masses(idx));
 cg = sum(([masses(idx);masses(idx)].*X(:,idx)),2)/tm;
 plot(cg(1),cg(2),'r*','markersize',7); 