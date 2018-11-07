set(groot,'defaultTextInterpreter','latex');
set(0, 'DefaultFigureVisible', 'on');
clc;close all;clear;
total_samples = [10000];
total_trials = 50;
per_list = [0:.19:1];
ref_param = [1,2];
n_measure_list = [2,3];
dist_name = 'lognormal';

n_str = sprintfc('%d',n_measure_list);
n_str = sprintfc('BMM %d pivots',2.^n_measure_list);

make_dist_handle = @(params)makedist(dist_name,params(1),params(2));
loss_handle_list ={};
for i_n = 1:length(n_measure_list)
    loss_handle_list{i_n} = @(x,per,params,tr)categorical_loss(make_dist_handle,x,per,params,n_measure_list(i_n));
end
% loss_handle_list{length(loss_handle_list)+1} =@(x,per,params,tr)TMLE_loss(x,params,per,make_dist_handle,-1,tr)
loss_handle_list{length(loss_handle_list)+1} =@(x,per,params,tr)TMLE_loss(x,params,per,make_dist_handle,0,tr)
loss_handle_list{length(loss_handle_list)+1} =@(x,per,params,tr)TMLE_loss(x,params,per,make_dist_handle,10,tr)
% loss_handle_list{length(loss_handle_list)+1} =@(x,per,params,tr)bmm(x,params,per,make_dist_handle,.5)
% loss_handle_list{length(loss_handle_list)+1} =@(x,per,params,tr)bmm(x,params,per,make_dist_handle,1,3)
% loss_handle_list{i_nc+1} = @(x,per,params,tr)loss(x,params(1),params(2),tr,per,0.1);
% loss_handle_list{i_n+2} = @(x,per,params,tr)loss(x,params(1),params(2),tr,per,0);
% n_str{length(n_str)+1} = 'TMLE REG';
n_str{length(n_str)+1} = 'CMLE';

n_str{length(n_str)+1} = 'CMAP';
% n_str{length(n_str)+1} = 'Binom Ordered';
% n_str{length(n_str)+1} = 'Binom';
% n_str{i_n+1} = 'TMLE';
% loss_handle_list = {@(x,per,params)categorical_loss(make_dist_handle,x,per,params,3)};
% per_list = 0:.3:1;

[err_mat,loss_mat] = test_measures(loss_handle_list,ref_param,make_dist_handle,per_list,total_trials,total_samples);

% wrong_guesses = sum(guess_mat,4)
% for i = 1:size(wrong_guesses,1)
%     guess_mat(i,:)./total_trials
%     figure;bar(categorical(n_str),wrong_guesses(i,:));
% end
err_mat_mean = mean(err_mat,4,'omitnan');
err_mat_std = std(err_mat,0,4,'omitnan');
loss_mat_mean = mean(loss_mat,4,'omitnan');
loss_mat_mean_per = squeeze(mean(loss_mat_mean,1,'omitnan'));
err_mat_mean_x_per = squeeze(mean(err_mat_mean,1,'omitnan'));
err_mat_std_x_per = squeeze(std(err_mat_mean,1,'omitnan'));
err_mat_mean_x_total_samples = squeeze(median(err_mat_mean,2,'omitnan'));

figure;plot(per_list,err_mat_mean_x_per');legend(n_str);grid on;xlabel('PER');ylabel(' $E_\theta := ||\hat{\theta}-\theta||$');title(sprintf('Full Distribution Total Samples %d, total trials: %d',total_samples,total_trials));saveas(gcf,sprintf('./Plots/Paper Plots/BME_Res%d %d %s.png',total_samples(1),total_trials,dist_name));
% figure;plot(per_list,squeeze(err_mat_std)');legend(n_str);grid on;xlabel('PER');ylabel(' $E_\theta := ||\hat{\theta}-\theta||$');title(sprintf('Full Distribution Total Samples %d, total trials: %d',total_samples,total_trials))
figure;plot(per_list,-loss_mat_mean_per');legend(n_str,'location','best');grid on;xlabel('PER');ylabel('Average Log Likelihood');title(sprintf('Average of loglikelihood (Training Total Samples %d, Total Trials: %d)',total_samples,total_trials));saveas(gcf,sprintf('./Plots/Paper Plots/BME_Res_LL%d %d %s.png',total_samples(1),total_trials,dist_name));
figure;
bin_edges = 0:1:50;
% for i = 1:size(loss_mat,3)
%     hold on;
%     histogram(squeeze(loss_mat(1,1,i,:)),bin_edges,'Normalization','cdf')
%     
% end
% legend(n_str)
% figure;
% bin_edges = 0:1:50;
% for i = 1:size(loss_mat,3)
%     figure;
%     histogram(squeeze(err_mat(1,1,i,:)),'Normalization','cdf')
%     
% end
% legend(n_str)
% saveas(gcf,'./Plots/Paper Plots/test_measures.png');
% plot(log10(total_samples),err_mat_mean_x_total_samples');legend(n_str);grid on;xlabel('$log(total samples)$');ylabel(' $E_\theta := ||\hat{\theta}-\theta||$');title(sprintf('Full Distribution Total Samples %d, total trials: %d',total_samples,total_trials))
