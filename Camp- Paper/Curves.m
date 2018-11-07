clc
close all
clear
set(groot,'defaultTextInterpreter','latex');
dist = makedist('weibull',3,2);
alpha = @(p)log(-log(1-p));
k = @(p,t,lambda)(alpha(p))./(log(t)-log(lambda));
p_array = [.1:.1:.9];
t = dist.icdf(p_array);
legend_str = sprintfc('F(%d)= %g t=%g',[1:length(t);p_array;t]');
lam = .6:.01:10;
for i = 1:length(p_array)
p_1 = p_array(i);
t_1 = t(i);
k_val_1 = k(p_1,t_1,lam);
lam_2 = lam(k_val_1>0);
k_val_1 = k_val_1(k_val_1>0);
plot(lam_2,k_val_1)
hold on
end

legend(legend_str,'location','best');
% lam = .6:.01:10;
% p_2 = .8;
% t_2 = dist.icdf(p_2);
% k_val_2 = k(p_2,t_2,lam);
% lam = lam(k_val_2>0);
% k_val_2 = k_val_2(k_val_2>0);
% plot(lam,k_val_2);
title('Weibull BM Function')
xlabel('$\lambda$');
ylabel('$K$');
grid on;
xlim([0.5,5])
ylim([0.5,5])
saveas(gcf,sprintf('./Plots/Paper Plots/Weibull_Bernouli_Curve.png'));