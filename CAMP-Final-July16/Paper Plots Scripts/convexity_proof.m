clc
close all
clear

%% Nakagami
mu = eps:.01:10;
omega = eps:.01:10;
fun_handle = @(x)log(pdf(makedist('normal',x,2),5));
mu_param = arrayfun(fun_handle,mu);
% fun_handle = @(x,y)pdf(makedist('nakagami',x,1),1);
figure;plot(mu,mu_param);