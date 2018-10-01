function [p_l,p_h] = categorical(cdf_handle,t)
%CATEGORICAL Summary of this function goes here
%   Detailed explanation goes here
p_h = 1-cdf_handle(t);
p_l = 1-p_h;
end

