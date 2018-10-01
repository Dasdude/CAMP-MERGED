function [loss] = categorical_loss(make_dist_handle,x,per,params,n)
%CATEGORICAL_LOSS Summary of this function goes here
%   Detailed explanation goes here
f = make_dist_handle(params);
cdf_handle = @(y)f.cdf(y);
max_x = max(x);
min_x =min(x);
% s_vals = ((max_x-min_x)*rand(1,n))+min_x;
total_samples = min(n,unique(length(x)));

s_vals = randsample(unique(x),total_samples);
s_vals = ((max_x-min_x)*([1:n-1]/(n)))+min_x;
% p_map = [1:n-1]
a = 1:n-1;
q = mod(a,[2.^(0:(log2(n)-1))]');
p_map = 2.^(-sum(q>0));
p_map = ones(size(s_vals));
% s_vals = [.125,.25,.375,.5,.625,.75,.875];
% p_map = [[1:n./2]./n,[(n/2) -1:-1:1]./n]
% s_vals = ((max_x-min_x)*s_vals)+min_x;
% p_map = [.125,.25,.125,.5,.125,.25,.125];
% s = sort(unique(x));
% prctile_array = 100*([1:n]./(n+1));
% s_vals = unique(prctile(x,prctile_array));
% points = min(n,length(s));
% if points<1
%     loss = 1e5;
%     return
% end
% i_0 = 1;
% piv = 1:points;
% i_list = int64(piv.*length(s)./(points));
% i_2 = int64(2*length(s)./points);
% index_list = i_list;
p_l = categorical_emperical(x,s_vals');
p_l_ref = p_l*(1-per) + per;
p_l_target = cdf_handle(s_vals);
loss = sum(p_map.*(-p_l_ref.*log(p_l_target)-(1-p_l_ref).*log(1-p_l_target)));
% for i = 1:length(s_vals)
%     p_l = categorical_emperical(x,s_vals(i));
%     p_l_ref = p_l*(1-per) + per;
%     p_l_target = cdf_handle(s_vals(i));
%     loss= loss-p_l_ref*log(p_l_target)-(1-p_l_ref)*log(1-p_l_target);
% end

end
