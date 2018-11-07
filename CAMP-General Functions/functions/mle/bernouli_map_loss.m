function [loss] = bernouli_map_loss(make_dist_handle,x,per,params,n)
%CATEGORICAL_LOSS Summary of this function goes here
%   Detailed explanation goes here
    if length(x)<4
        loss=0;
        return
    end
    f = make_dist_handle(params);
    cdf_handle = @(y)f.cdf(y);
    max_x = max(x);
    min_x =min(x);
    s_vals = ((max_x-min_x)*([0:n-1]/(n)))+min_x;
    s_vals = unique(x);
    p_l = categorical_emperical(x,s_vals');
    p_l_ref = p_l*(1-per) + per;
    p_l_target = max(eps,min(1-eps,cdf_handle(s_vals)));
    h_ref = -(p_l_ref.*log2(min(max(p_l_ref,eps),1-eps))-(1-p_l_ref).*log2(min(max(1-p_l_ref,eps),1-eps)));
    loss = mean((2.^(h_ref)).*(-p_l_ref.*log2(p_l_target)-(1-p_l_ref).*log2(1-p_l_target)));
end

