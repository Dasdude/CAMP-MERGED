function [loss] = median_qnt_loss_emp(x,p_min,p_max,n,f_handle,max_n,lambda)
%MEDIAN_QUANTILE_LOSS Summary of this function goes here
%   Detailed explanation goes here

length(x);
if length(x)==1||length(unique(x))==1||n==max_n+1
    loss=0;
    return
end
t = median(x);
x_l = x(x<t);
x_h = x(x>t);
x = x(x~=t);
sprintf('%d','%d',length(x_l),length(x_h));
%%PRIOR HERE
p_emp = (length(x_l))./(length(x));
p_emp = min(p_emp,1-eps);
p_emp = max(eps,p_emp);
p_mid = f_handle.cdf(t);
p_mid = min(p_mid,1-eps);
p_l = (p_mid-p_min)./(p_max-p_min+eps);
p_h =1-p_l;
loss = -(p_emp*(log2(p_l))+((1-p_emp)*(log2(p_h))));
% h_emp = p_emp*(log2(p_emp))+((1-p_emp)*(log2(1-p_emp)));
% h_f_emp = p_l*(log2(p_l))+((p_h)*(log2(p_h)));
if  isnan(loss)||p_mid<p_min
%     sprintf('stopped at n: %d',n)
    loss=0;
    return;
end

h_l = median_qnt_loss_emp(x_l,p_min,p_mid,n+1,f_handle,max_n,lambda);
h_h = median_qnt_loss_emp(x_h,p_mid,p_max,n+1,f_handle,max_n,lambda);
loss = loss+lambda*(p_emp*(h_l)+(1-p_emp)*(h_h));

