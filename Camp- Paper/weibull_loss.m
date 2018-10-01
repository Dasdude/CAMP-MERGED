clc
clear
close all

load data_sample.mat
lam = .3:.05:3;
k = .2:.05:2.5;
x = fading_linear_cell{100};
% x =rand(1,20);
[lam_mesh,k_mesh,x_mesh] = meshgrid(lam,k,x);
loss = @(x,m,k)-(log(k)- k.*log(m) +(k-1).*log(x)-(x./m).^k);
loss_map = sum(loss(x_mesh,lam_mesh,k_mesh),3);
[lk_mesh,kl_mesh] = meshgrid(lam,k);
loss_map(loss_map>1e4)=nan;
figure;surf(lk_mesh,kl_mesh,loss_map);xlabel('\lambda');ylabel('k');zlabel('loss');xlim([min(lam),max(lam)]);ylim([min(k),max(k)]);
% hold on; plot(

[xk_mesh,kx_mesh] = meshgrid(x,k);
figure;opt_lambda = (sum(xk_mesh.^(kx_mesh-1),2)).^(1./k');plot(k,opt_lambda);

% dllam = @(lam,k)