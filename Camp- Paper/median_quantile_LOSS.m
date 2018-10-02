function [loss] = median_quantile_LOSS(x,p,n,f_handle,max_n)
%MEDIAN_QUANTILE_LOSS Summary of this function goes here
%   Detailed explanation goes here
length(x);
if length(x)==1||length(unique(x))==1||n==max_n
    loss=0;
    return
end
t = f_handle.icdf(p);
x_l = x(x<t);
x_h = x(x>t);
sprintf('%d','%d',length(x_l),length(x_h));
p_emp = (length(x_l))./(length(x));
p_emp = min(p_emp,1-eps);
p_emp = max(eps,p_emp);
loss = 1+ p_emp.*log2(p_emp)+(1-p_emp).*log2(1-p_emp);
if  isnan(loss)
    sprintf('stopped at n: %d',n)
    loss=0;
    return;
end
% f_handle_l = @(x)(f_handle(.5- (x/2)) );
% f_handle_h = @(x)(icdf_handle(.5+ (x/2)) );
h_l = median_quantile_LOSS(x_l,p-(.5).^(n+1),n+1,f_handle,max_n);
h_h = median_quantile_LOSS(x_h,p+(.5).^(n+1),n+1,f_handle,max_n);
loss = loss+.5*(p_emp*h_l+(1-p_emp)*h_h);

