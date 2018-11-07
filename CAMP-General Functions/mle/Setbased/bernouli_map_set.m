function [loss] = bernouli_map_set(make_dist_handle,x,per,params,n)
%CATEGORICAL_LOSS Summary of this function goes here
%   Detailed explanation goes here
ll_list = zeros(1,length(x));
samples_per_cell = funoncellarray1input(x,@length);
samples_per_cell(isnan(samples_per_cell))=0;
for i = 1:length(x)
    if isempty(x{i})
        per(i)=0;
        continue;
    end
    ll_list(i)=bernouli_map_loss(make_dist_handle,x{i}',per(i),params,n);
    
end
loss = mean(ll_list.*(samples_per_cell./(1-per'))./sum((samples_per_cell./(1-per'))));
end

