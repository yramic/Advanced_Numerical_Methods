figure('Name','Gaussweights');
hold on;
[x,w] = fclencurt(2,-1,1); plot(x,w,'r+'); sum(w)
[x,w] = fclencurt(4,-1,1); plot(x,w,'bo');
[x,w] = fclencurt(6,-1,1); plot(x,w,'c*');
[x,w] = fclencurt(8,-1,1); plot(x,w,'mx');
[x,w] = fclencurt(10,-1,1); plot(x,w,'k^');
[x,w] = fclencurt(12,-1,1); plot(x,w,'gp');
[x,w] = fclencurt(14,-1,1); plot(x,w,'rh'); sum(w)
set(gca,'fontsize',12);
title('{\bf Clenshaw-Curtis weights for [-1,1]}','FontSize',20);
xlabel('{t_j}','fontsize',14);
ylabel('{\bf w_j}','fontsize',14);
legend('n=2','n=4','n=6','n=8','n=10','n=12','n=14');
axis([-1 1 0 1.1]);
hold off;

print -depsc2 'CCweights.eps';
