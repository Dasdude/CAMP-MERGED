function [p] = sigmoid(density,a,x_bias)
%SIGMOID Give Density and params provided and the output is the weight used
%for high density. the parmaeters for the density 
%   Detailed explanation goes here
% p = p_fun(density,a,x_bias)
    f = @(x,a)exp(a*(x-x_bias))./(1+exp(a*(x-x_bias)));
    if isempty(density)
        x=1:60;
        figure;plot(x,f(x,a));
        return
    end
    p = f(density,a);
    
%     figure;plot(x,f(x,a));
end

