figure('Name','CC nodes');
hold on;
for i=2:20
  plot(fclencurt(i,-1,1),i*ones(1,i),'rp');
end
set(gca,'fontsize',12);
title('Clenshaw-Curtis nodes in [-1,1]','FontSize',20);
xlabel('t');
ylabel('{\bf Number n of quadrature nodes}','fontsize',14);
axis([-1 1 1.5 20.5]);
hold off;

print -depsc2 'CCnodes.eps';