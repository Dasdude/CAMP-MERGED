% clc
% close all
% clear

f_handle=@(theta) makedist('weibull',theta(1),theta(2));
x = random(f_handle([1,1]),1,1e1);
% loss_median = median_quantile_LOSS(x,.5,2,f_handle,10);
% f_handle= @(theta)makedist('weibull',theta(1),theta(2));
loss_ll = TMLE_loss(x,[1,1],0,f_handle,0,0);

f_median =@(x,per,params,tr)bmm(x,params,per,f_handle,.5)
% f_median(x,0,[1,2],0)
% f_median(x,0,[1,1],0)
% f_median(x,0,[1,1.1],0)
fun = @(theta)(f_median(x,0,theta,0));
% fun_median = @(theta)median_quantile_LOSS(x,.5,2,f_handle,10);
options = optimoptions(@fmincon,'PlotFcns',@optimplotfval,'Display','iter','Algorithm','interior-point','MaxFunctionEvaluations',10000);
% options = optimset('PlotFcns',@optimplotfval);
[params,loss_val,~,out]=fmincon(fun,[1,2],[],[],[],[],[eps,eps],[],[],options);
% [params,loss_val,~,out]=fminsearch(fun,[1,2],options);
sol_dist = makedist('weibull',params(1),params(2));
a = histogram(sol_dist.random(1,1e4),'normalization','cdf');hold on;histogram(x,'normalization','cdf','BinEdges',a.BinEdges);
params