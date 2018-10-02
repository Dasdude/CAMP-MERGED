function [loss] = bmm(x,theta,per,f_handle,lambda,max_n)
%BMM Summary of this function goes here
%   Detailed explanation goes here
ft_handle = f_handle(theta);
% x,p_min,p_max,n,f_handle,max_n
% loss = median_quantile_LOSS(x,.5,1,ft_handle,max_n);
max_p = 1;
loss = median_qnt_loss_emp(x,0,max_p,1,ft_handle,max_n,lambda);
end

