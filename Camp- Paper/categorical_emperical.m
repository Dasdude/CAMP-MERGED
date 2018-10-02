function [p_l] = categorical_emperical(x,t)
%CATEGORICAL_EMPERICAL Summary of this function goes here
%   Detailed explanation goes here
flag = x>t;
p_h = sum(flag')./length(x);
p_l = 1-p_h;
end

