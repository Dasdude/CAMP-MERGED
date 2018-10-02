set(groot,'defaultTextInterpreter','latex');
clc;close all;clear;
total_samples = [1e1,50,100,1000,5000];
total_trials = 300;
per_list = [0:.05:1];
ref_param = [.8,1.2];
n_measure_list = [4,8,32,64,128,512];
f = @(x,m,k) ((k./m).*((x./m).^(k-1)).*exp(-((x./m).^k)));
f_cdf = @(x,m,k) 1-exp(-(x./m).^k);
icdf_f = @(p,m,k) m.*(-log(1-p)).^(1./k);
p_hat = @(x,m,k,tr)(f(x,m,k)./(1-f_cdf(tr,m,k))).*((x>=tr));
reg = @(x,m,k,tr,p)(icdf_f(p,m,k)-tr).^2;
loss = @(x,m,k,tr,p,lambda)-squeeze(mean(log(p_hat(x,m,k,tr)),2))+lambda.*squeeze(mean(reg(x,m,k,tr,p),2));

dist_name = 'nakagami';

n_str = sprintfc('%d',n_measure_list);

make_dist_handle = @(params)makedist(dist_name,params(1),params(2));
loss_handle_list ={};
for i_n = 1:length(n_measure_list)
    loss_handle_list{i_n} = @(x,per,params)categorical_loss(make_dist_handle,x,per,params,n_measure_list(i_n));
end
loss_handle_list{i_n+1} =@(x,per,params,tr)TMLE_loss(x,params,per,make_dist_handle,.1,tr)
loss_handle_list{i_n+2} =@(x,per,params,tr)TMLE_loss(x,params,per,make_dist_handle,0,tr)
% loss_handle_list{i_n+1} = @(x,per,params,tr)loss(x,params(1),params(2),tr,per,0.1);
% loss_handle_list{i_n+2} = @(x,per,params,tr)loss(x,params(1),params(2),tr,per,0);
n_str{i_n+1} = 'TMLE REG';
n_str{i_n+2} = 'TMLE';
% n_str{i_n+1} = 'TMLE';
% loss_handle_list = {@(x,per,params)categorical_loss(make_dist_handle,x,per,params,3)};
% per_list = 0:.3:1;

[err_mat,param_mat] = test_measures(loss_handle_list,ref_param,make_dist_handle,per_list,total_trials,total_samples);
err_mat_mean = squeeze(median(err_mat,4,'omitnan'));
err_mat_mean_x_per = squeeze(median(err_mat_mean,1,'omitnan'));
err_mat_mean_x_total_samples = squeeze(median(err_mat_mean,3,'omitnan'));
plot(per_list,err_mat_mean_x_per');legend(n_str);grid on;xlabel('PER');ylabel(' $E_\theta := ||\hat{\theta}-\theta||$');title(sprintf('Full Distribution Total Samples %d, total trials: %d',total_samples,total_trials))
plot(log10(total_samples),err_mat_mean_x_total_samples');legend(n_str);grid on;xlabel('$log(total samples)$');ylabel(' $E_\theta := ||\hat{\theta}-\theta||$');title(sprintf('Full Distribution Total Samples %d, total trials: %d',total_samples,total_trials))
