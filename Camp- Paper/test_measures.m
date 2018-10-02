function [err_mat,param_mat] = test_measures(loss_handle_list,param_ref,make_dist_handle,per_samples_list,total_trials,total_samples_list)
%TEST_MEASURES Summary of this function goes here
%   Detailed explanation goes here
f_ref = make_dist_handle(param_ref);
err_mat  = zeros(length(total_samples_list),length(per_samples_list),length(loss_handle_list),total_trials)*nan;
param_mat = zeros(length(total_samples_list),length(per_samples_list),length(loss_handle_list),total_trials,length(param_ref))*nan;
log_like = zeros(length(total_samples_list),length(per_samples_list),length(loss_handle_list),total_trials,length(param_ref))*nan;
% total_samples = total_samples_list;
for total_samples_idx = 1:length(total_samples_list)
    total_samples = total_samples_list(total_samples_idx);
    for per_samples_idx = 1:length(per_samples_list)
        per_samples = per_samples_list(per_samples_idx);
        for i = 1:total_trials
            x = random(f_ref,[total_samples,1]);
    %         x = dbm2linear(double(int64(linear2dbm(x))));
            for loss_idx = 1:length(loss_handle_list)
                clc
                sprintf('Total_samples: %d PER: %d trial:%d/%d loss_idx:%d/%d ',total_samples,per_samples,i,total_trials,loss_idx,length(loss_handle_list))
                tr = icdf(f_ref,per_samples);
                if loss_idx >= length(loss_handle_list)-1
                    loss_fun = loss_handle_list{loss_idx};
                    loss_fun = @(a,b,c)loss_fun(a,b,c,tr);
                else
                    loss_fun = loss_handle_list{loss_idx};
                end
                x_trunc= x(x>tr);
                if length(x_trunc)<5
                    continue
                end
                per = 1-(length(x_trunc)./length(x));

                if per==1
                    continue
                end
                fun = @(theta)loss_fun(x_trunc',per,theta);
                options = optimoptions(@fmincon,'Display','off','Algorithm','interior-point','MaxFunctionEvaluations',10000);
                [params,loss_val,~,out]=fmincon(fun,[1,1],[],[],[],[],[.5,eps],[],[],options);
                if any(isnan(params))
                    disp('something is wrong'); 
                end
                err_mat(total_samples_idx,per_samples_idx,loss_idx,i) =sqrt(sum((params-param_ref).^2));
                param_mat(total_samples_idx,per_samples_idx,loss_idx,i,:) =params;
             end
        end
    end
end
end

