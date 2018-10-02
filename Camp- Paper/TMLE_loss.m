function [loss] = TMLE_loss(x,theta,per,f_handle,lambda,t)
%TMLE_LOSS Summary of this function goes here
%   Detailed explanation goes here
f_handle = f_handle(theta);
p_hat = f_handle.pdf(x)./(1-f_handle.cdf(t)+eps);
p = f_handle.pdf(x);
if lambda<0
    p_emp = categorical_emperical(x,t); 
%     loss = -((mean(log(p_hat+eps)))+per*(p_emp*log(f_handle.cdf(t)+eps)+(1-p_emp)*log(1-f_handle.cdf(t)+eps)));
    loss = -((1-per)*(sum(log(p_hat+eps))+per*(p_emp*log(f_handle.cdf(t)+eps))));
else
    loss = -mean(log(p_hat+eps))+ lambda.*abs(t-f_handle.icdf(per+eps));
end
end

