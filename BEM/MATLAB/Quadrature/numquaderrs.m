function numquaderrs()
% Numerical quadrature on [0,1]
% Numerical experiments for Sect 4.1 of ADVNCSE lecture material
N = 20;

figure('Name','tent');
exact = 0.25;
chbres = numquad(@tentfn,0.0,1,N,'CC');
gaures = numquad(@tentfn,0.0,1,N,'Gauss');
loglog(chbres(:,1),abs(chbres(:,2)-exact),'m+-',...
       gaures(:,1),abs(gaures(:,2)-exact),'r+-');
set(gca,'fontsize',12);
title('Numerical quadrature of tent function');
xlabel('{\bf Number of quadrature nodes}','fontsize',14);
ylabel('{\bf |quadrature error|}','fontsize',14);
legend('Clenshaw-Curtis quadrature',...
       'Gauss quadrature','location','best');

print -depsc2 'globquad_tent.eps';

figure('Name','bump');
princ = @(x) 0.5*x - 0.0625*sin(4 - 8*x);
exact = princ(0.75)-princ(0.25);
chbres = numquad(@cosbump,0.0,1,N,'CC');
gaures = numquad(@cosbump,0.0,1,N,'Gauss');
loglog(chbres(:,1),abs(chbres(:,2)-exact),'m+-',...
       gaures(:,1),abs(gaures(:,2)-exact),'r+-');
set(gca,'fontsize',12);
title('Numerical quadrature of bump function');
xlabel('{\bf Number of quadrature nodes}','fontsize',14);
ylabel('{\bf |quadrature error|}','fontsize',14);
legend('Clenshaw-Curtis quadrature',...
       'Gauss quadrature','location','best');

print -depsc2 'globquad_bump.eps';

figure('Name','t*log t');
princ = @(x) 0.25*x^2*(2*log(x)-1);
exact = princ(1.0);
eqdres = numquad(@xlogx,0.0,1,N,'equidistant');
chbres = numquad(@xlogx,0.0,1,N,'CC');
gaures = numquad(@xlogx,0.0,1,N,'Gauss');
loglog(eqdres(:,1),abs(eqdres(:,2)-exact),'b+-',...
       chbres(:,1),abs(chbres(:,2)-exact),'m+-',...
       gaures(:,1),abs(gaures(:,2)-exact),'r+-');
set(gca,'fontsize',12);
title('Numerical quadrature of function  t*log t');
xlabel('{\bf Number of quadrature nodes}','fontsize',14);
ylabel('{\bf |quadrature error|}','fontsize',14);
legend('Equidistant Newton-Cotes quadrature',...
       'Clenshaw-Curtis quadrature',...
       'Gauss quadrature','location','best');
eqrerr = abs(eqdres(:,2)-exact);
chberr = abs(chbres(:,2)-exact);
gauerr = abs(gaures(:,2)-exact);
eqdp1 = polyfit(eqdres(:,1),log(eqrerr),1)
chbp1 = polyfit(chbres(:,1),log(chberr),1)
gaup1 = polyfit(gaures(:,1),log(gauerr),1)
grid on;

print -depsc2 'globquad_tlogt.eps';

figure('Name','log(t+0.1)');
logprinc = @(x) x*(log(x)-1);
exact = logprinc(1.1)-logprinc(0.1);
eqdres = numquad(inline('log(x+0.1)'),0,1,N,'equidistant');
chbres = numquad(inline('log(x+0.1)'),0,1,N,'CC');
gaures = numquad(inline('log(x+0.1)'),0,1,N,'Gauss');
semilogy(eqdres(:,1),abs(eqdres(:,2)-exact),'b+-',...
     chbres(:,1),abs(chbres(:,2)-exact),'m+-',...
     gaures(:,1),abs(gaures(:,2)-exact),'r+-');
set(gca,'fontsize',12);
title('Numerical quadrature of function  log(t+0.1)');
xlabel('{\bf Number of quadrature nodes}','fontsize',14);
ylabel('{\bf |quadrature error|}','fontsize',14);
legend('Equidistant Newton-Cotes quadrature',...
       'Clenshaw-Curtis quadrature',...
       'Gauss quadrature','location','best');
eqdp1 = polyfit(eqdres(:,1),log(abs(eqdres(:,2)-exact)),1)
chbp1 = polyfit(chbres(:,1),log(abs(chbres(:,2)-exact)),1)
gaup1 = polyfit(gaures(:,1),log(abs(gaures(:,2)-exact)),1)
grid on;

print -depsc2 'globquad_logt.eps';

figure('Name','1/(1+(5t)^2)');
exact = atan(5)/5;
eqdres = numquad(inline('1./(1+(5*x).^2)'),0,1,N,'equidistant');
chbres = numquad(inline('1./(1+(5*x).^2)'),0,1,N,'CC');
gaures = numquad(inline('1./(1+(5*x).^2)'),0,1,N,'Gauss');
semilogy(eqdres(:,1),abs(eqdres(:,2)-exact),'b+-',...
     chbres(:,1),abs(chbres(:,2)-exact),'m+-',...
     gaures(:,1),abs(gaures(:,2)-exact),'r+-');
set(gca,'fontsize',12);
title('Numerical quadrature of function  log(t+0.1)');
xlabel('{\bf Number of quadrature nodes}','fontsize',14);
ylabel('{\bf |quadrature error|}','fontsize',14);
legend('Equidistant Newton-Cotes quadrature',...
       'Clenshaw-Curtis quadrature',...
       'Gauss quadrature','location','best');
eqdp1 = polyfit(eqdres(:,1),log(abs(eqdres(:,2)-exact)),1)
chbp1 = polyfit(chbres(:,1),log(abs(chbres(:,2)-exact)),1)
gaup1 = polyfit(gaures(:,1),log(abs(gaures(:,2)-exact)),1)
grid on;

print -depsc2 'globquad_runge.eps';

figure('Name','sqrt(t)');
exact = 2/3;
eqdres = numquad(inline('sqrt(x)'),0,1,N,'equidistant');
chbres = numquad(inline('sqrt(x)'),0,1,N,'CC');
gaures = numquad(inline('sqrt(x)'),0,1,N,'Gauss');
loglog(eqdres(:,1),abs(eqdres(:,2)-exact),'b+-',...
     chbres(:,1),abs(chbres(:,2)-exact),'m+-',...
     gaures(:,1),abs(gaures(:,2)-exact),'r+-');
set(gca,'fontsize',12);
axis([1 25 0.000001 1]);
title('Numerical quadrature of function sqrt(t)');    
xlabel('{\bf Number of quadrature nodes}','fontsize',14);
ylabel('{\bf |quadrature error|}','fontsize',14);
legend('Equidistant Newton-Cotes quadrature',...
       'Clenshaw-Curtis quadrature',...
       'Gauss quadrature','location','best');
eqdp2 = polyfit(log(eqdres(:,1)),log(abs(eqdres(:,2)-exact)),1)
chbp2 = polyfit(log(chbres(:,1)),log(abs(chbres(:,2)-exact)),1)
gaup2 = polyfit(log(gaures(:,1)),log(abs(gaures(:,2)-exact)),1)
grid on;

print -depsc2 'globquad_root.eps';



