function plotpartitionlevel(X,masses,l)
% Plot stars and clustering boxes at a certain level

d = 2^(-l); % Size of boxes
figure('name','Clustering boxes');
plot([0 1 1 0 0],[0 0 1 1 0],'k-','linewidth',2); hold on;
title(sprintf('Cluster boxes at level %i',l));
plot(X(1,:),X(2,:),'b.');
for j=1:2^l
    for k=1:2^l
        clusterplot([(j-1)*d,j*d,(k-1)*d,k*d],X,masses);
    end
end
hold off;