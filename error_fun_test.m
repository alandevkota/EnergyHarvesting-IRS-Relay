x=0:0.01:2;
f1=erf(x);
f2=1-exp(-x.^2)./(sqrt(pi).*x);

figure
plot(x,f1,'--','color',[1 0 0],'LineWidth',2,'MarkerSize',9);hold on;
plot(x,f2,'-','color',[0 0 1],'LineWidth',2);hold on;
ylim([0 2])
