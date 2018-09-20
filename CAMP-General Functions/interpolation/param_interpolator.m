function [param_results] = param_interpolator(density,low_density_param,high_density_param)
%PARAM_INTERPOLATOR Summary of this function goes here
%   Detailed explanation goes here
a= .3;
x_bias = 22;
f = @(x,a)exp(a*(x-x_bias))./(1+exp(a*(x-x_bias)));
p = f(density,a);
param_results = p*high_density_param+((1-p)*(low_density_param));
end

