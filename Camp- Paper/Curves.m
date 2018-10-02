clc
close all
clear
set(groot,'defaultTextInterpreter','latex');
dist = makedist('weibull',3,2);
alpha = @(p)log(-log(1-p));
k = @(p,t,lambda)(alpha(p))./(log(t)-log(lambda));
lam = .6:.01:10;
p_1 = .5;
t_1 = dist.icdf(p_1);
k_val_1 = k(p_1,t_1,lam);
lam = lam(k_val_1>0);
k_val_1 = k_val_1(k_val_1>0);
plot(lam,k_val_1)
hold on
lam = .6:.01:10;
p_2 = .8;
t_2 = dist.icdf(p_2);
k_val_2 = k(p_2,t_2,lam);
lam = lam(k_val_2>0);
k_val_2 = k_val_2(k_val_2>0);
plot(lam,k_val_2);legend(sprintf('F(t_1)= %g t_1=%g',p_1,t_1),sprintf('F(t_2)= %g t_2=%g',p_2,t_2),'location','best');
title('Weibull BM Function')
xlabel('$\lambda$');
ylabel('$K$');
grid on;
xlim([0.5,5])
ylim([0.5,5])
saveas(gcf,sprintf('./Plots/Paper Plots/Weibull_Bernouli_Curve.png'));