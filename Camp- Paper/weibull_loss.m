clc
clear
close all
addpath(genpath('.'))
load data_sample.mat
lam = .3:.05:5;
k = .2:.05:4;
per = 0.9;
fading_linear_cell = dbm2linear(dataset_cell);
ref_param = [3,2.5];
n_measure_list =2.^[1:8];
total_samples = 50

make_dist_handle = @(x)makedist('weibull',x(1),x(2));
ref_dist = make_dist_handle(ref_param);
x = random(make_dist_handle(ref_param),1,total_samples);
% x = dbm2linear(floor(linear2dbm(x)));
t = icdf(ref_dist,per);
x_trunc =x(x>=t);
per = 1- (length(x_trunc)./length(x));
n_str = sprintfc('%d',n_measure_list);
loss_handle_list = {};
for i_n = 1:length(n_measure_list)
    loss_handle_list{i_n} = @(x,m,k)categorical_loss(make_dist_handle,x,per,[m,k],n_measure_list(i_n));
end
loss_handle_list{length(loss_handle_list)+1} =@(x,m,k)TMLE_loss(x,[m,k],per,make_dist_handle,0,t);n_str{length(n_str)+1} = 'TMLE';

for i = 1:length(loss_handle_list)
% loss= @(x,m,k)categorical_loss(makedist_handle,x_trunc,per,[m,k],8);
loss_val = zeros(length(lam),length(k));
loss = loss_handle_list{i};

for lam_idx = 1:length(lam)
    for k_idx = 1:length(k)
        m = lam(lam_idx);
        k_val = k(k_idx);
        loss_val(lam_idx,k_idx) = loss(x_trunc,m,k_val);
    end
end
loss_val = loss_val';
% loss_val(isnan(loss_val))=inf;
% loss_val(loss_val>.7)=nan;

[lk_mesh,kl_mesh] = meshgrid(lam,k);
% figure;surf(lk_mesh,kl_mesh,loss_val);xlabel('\lambda');ylabel('k');zlabel('loss');xlim([min(lam),max(lam)]);ylim([min(k),max(k)]);
vals = prctile(loss_val(:),2.^[-1:5]);
lam_vec = lk_mesh(:);
k_vec  = kl_mesh(:);
lam_sol = lam_vec(loss_val(:)==min(loss_val(:)));
k_sol = k_vec(loss_val(:) == min(loss_val(:)));
% lam_vec(min(loss_val(:))
figure;contour(lk_mesh,kl_mesh,loss_val,vals);hold on;scatter(ref_param(1),ref_param(2));hold on;scatter(lam_sol,k_sol);
title(n_str{i})
colorbar()
end
