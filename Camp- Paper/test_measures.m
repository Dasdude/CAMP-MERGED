function [err_mat,loss_mat] = test_measures(loss_handle_list,param_ref,make_dist_handle,per_samples_list,total_trials,total_samples_x)
%TEST_MEASURES Summary of this function goes here
%   Detailed explanation goes here
f_ref = make_dist_handle(param_ref);
err_mat  = zeros(length(total_samples_x),length(per_samples_list),length(loss_handle_list),total_trials)*nan;
param_mat = zeros(length(total_samples_x),length(per_samples_list),length(loss_handle_list),total_trials,length(param_ref))*nan;
guess_mat = zeros(length(total_samples_x),length(per_samples_list),length(loss_handle_list),total_trials)*nan;
loss_mat = zeros(length(total_samples_x),length(per_samples_list),length(loss_handle_list),total_trials)*nan;
param_init = param_ref;
% param_init = [1,1]
% total_samples = total_samples_list;
samples_general = random(make_dist_handle(param_ref),[1e5,1]);
for total_samples_idx = 1:length(total_samples_x)
    total_samples = total_samples_x(total_samples_idx);
%     x_full = total_samples
    for per_samples_idx = 1:length(per_samples_list)
        per_samples = per_samples_list(per_samples_idx);
        for i = 1:total_trials
            
%             f_ref = make_dist_handle(rand(1,2)*5+1);
            
            x = random(f_ref,[total_samples,1]);
    %         x = dbm2linear(double(int64(linear2dbm(x))));
            tr = icdf(f_ref,per_samples);
            x_trunc= x(x>tr);
            per = 1-(length(x_trunc)./length(x));
            for loss_idx = 1:length(loss_handle_list)
                clc
                sprintf('Total_samples: %d PER: %d trial:%d/%d %d  ',total_samples,per_samples,i,total_trials,loss_idx)    
                loss_fun = loss_handle_list{loss_idx};
                loss_fun = @(a,b,c)loss_fun(a,b,c,tr);                
                fun = @(theta)loss_fun(x_trunc',per,theta);
                if per==1
                    continue
                end
                
                options = optimoptions(@fmincon,'Display','off','Algorithm','interior-point','MaxFunctionEvaluations',10000);
                [params,loss_val,~,out]=fmincon(fun,param_init,[],[],[],[],[eps,eps],[],[],options);
%                 params = [params,param_ref(2)];
                if any(isnan(params))
                    disp('something is wrong'); 
                end
                err_mat(total_samples_idx,per_samples_idx,loss_idx,i) =sqrt(sum((params-param_ref).^2));
                param_mat(total_samples_idx,per_samples_idx,loss_idx,i,:) =params;
                d_ref = make_dist_handle(params);
                loss_mat(total_samples_idx,per_samples_idx,loss_idx,i) = -mean(log2(d_ref.pdf(samples_general)+eps));
%                 loss_mat(total_samples_idx,per_samples_idx,loss_idx,i) =loss_val;
            end
%             param_dif =  squeeze(param_mat(total_samples_idx,per_samples_idx,:,i,:))-param_ref ;
%             param_err = sum(param_dif.^2,2);
%             [~,idx]=min(param_err);
%             param_truth = squeeze(param_mat(total_samples_idx,per_samples_idx,idx,i,:));
%             for loss_idx = 1:length(loss_handle_list)
%                 loss_fun = loss_handle_list{loss_idx};
%                 loss_fun = @(a,b,c)loss_fun(a,b,c,tr);                
%                 fun = @(theta)loss_fun(x_trunc',per,theta);
%                 loss_truth = fun(param_truth);
%                 guess_mat(total_samples_idx,per_samples_idx,loss_idx,i) = loss_mat(total_samples_idx,per_samples_idx,loss_idx,i) <loss_truth;
%             end
%             squeeze(guess_mat(total_samples_idx,per_samples_idx,:,i));
        end
    end
end
end

