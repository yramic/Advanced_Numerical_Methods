% Computation of axis of inertia
function axisofinertia()
% number of points
N = 100;
% Initialize points (random)
r = rand(2,N);
pts = [r(1,:).*cos(2*pi*r(2,:));2*r(1,:).*sin(2*pi*r(2,:))];
% center of gravity
c = sum(pts,2)/N;
% Inertia matrix
d = pts-repmat(c,1,N); M = d*d';
% diagonalize M, do sth. equivalent by means of SVD
[V,D] = eig(M);
[V,S,U] = svd(d,'econ');
% Columns of V give the axes of inertia
figure('name','Inertial splitting');
plot(pts(1,:),pts(2,:),'b+'); hold on;
a = c+V(:,1); b = c-V(:,1);
plot([a(1),b(1)],[a(2),b(2)],'r-');
a = c+V(:,2); b = c-V(:,2);
plot([a(1),b(1)],[a(2),b(2)],'m-');

pts1 = pts(:,find(((U(:,1) > 0) .* (U(:,2) > 0)) >0));
plot(pts1(1,:),pts1(2,:),'r^');
pts1 = pts(:,find(((U(:,1) <= 0) .* (U(:,2) > 0)) >0));
plot(pts1(1,:),pts1(2,:),'k^');
pts1 = pts(:,find(((U(:,1) > 0) .* (U(:,2) <= 0)) >0));
plot(pts1(1,:),pts1(2,:),'r*');
pts1 = pts(:,find(((U(:,1) <= 0) .* (U(:,2) <= 0)) >0));
plot(pts1(1,:),pts1(2,:),'k*');