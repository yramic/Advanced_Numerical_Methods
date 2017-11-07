function grav_clust_driver()

n = 200; % Number of start
% Distribute stars randomly in the unit square
X = rand(2,n);
% Slightly varying masses
masses = 1+rand(1,n);

% Top level cluster boxes
plotpartitionlevel(X,masses,1);
print -depsc2 'starclusterbox1.eps';
plotpartitionlevel(X,masses,2);
print -depsc2 'starclusterbox2.eps';
plotpartitionlevel(X,masses,3);
print -depsc2 'starclusterbox3.eps';

% Depict a special situation
figure('name','Admissible boxes');
plot([0 1 1 0 0],[0 0 1 1 0],'k-','linewidth',2); hold on;
plot(X(1,:),X(2,:),'c.');
title('Admissible clusters for a single star');
plot(0.9,0.15,'r.','markersize',12);
clusterplot([0 0.5 0.5 1],X,masses);
clusterplot([0 0.25 0 0.25],X,masses);
clusterplot([0 0.25 0.25 0.5],X,masses);
clusterplot([0.25 0.5 0 0.25],X,masses);
clusterplot([0.25 0.5 0.25 0.5],X,masses);
clusterplot([0.5 0.75 0.5 0.75],X,masses);
clusterplot([0.5 0.75 0.75 1],X,masses);
clusterplot([0.75 1 0.5 0.75],X,masses);
clusterplot([0.75 1 0.75 1],X,masses);
clusterplot([0.5 0.75 0.25 0.5],X,masses);
clusterplot([0.5 0.625 0 0.125],X,masses);
clusterplot([0.5 0.625 0.125 0.25],X,masses);
clusterplot([0.75 0.875 0.375 0.5],X,masses);
clusterplot([0.875 1 0.375 0.5],X,masses);
clusterplot([0.625 0.75 0.125 0.25],X,masses);
clusterplot([0.75 0.875 0.25 0.375],X,masses);
clusterplot([0.625 0.75 0 0.125],X,masses);
print -depsc2 'starclustadm.eps';

