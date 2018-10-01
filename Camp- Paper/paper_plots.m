clear
load data_sample.mat
close all
plot_folder = fullfile('Plots','Paper Plots');
mkdir(plot_folder)
data_fading_dbm = extract_fading(dataset_cell,pathloss,14);
fading_linear_cell = dbm2linear(data_fading_dbm);
generated_fading_linear = sample_generator(dist_obj,fading_params,1e3);
generated_fading_dbm = linear2dbm(generated_fading_linear);
generated_rssi_dbm = add_fading(pathloss,generated_fading_dbm,14);
[generated_rssi_dbm_truncated,generated_per,gen_pl_stat] = censor_data(generated_rssi_dbm,censor_function_handle);

%% Weibull Loss
% x = fading_linear_cell{200};
% t = dbm2linear(-90+pathloss(200)-14);
% lam =eps:.05:3;
% k = eps:.05:3;
% [lam_mesh,k_mesh,x_mesh] = meshgrid(lam,k,x);
% f = @(x,m,k) ((k./m).*((x./m).^(k-1)).*exp(-((x./m).^k)));
% f_cdf = @(x,m,k) 1-exp(-(x./m).^k);
% icdf_f = @(p,m,k) m.*(-log(1-p)).^(1./k);
% f_hat =@(x,m,k,tr) (f(x,m,k)./(1-f_cdf(tr,m,k))).*((x>=tr))
% loss_hat = @(x,m,k,tr)-log(f_hat(x,m,k,tr));
% loss = @(x,m,k)-log(f(x,m,k));
% % loss_map = mean(loss(x_mesh,lam_mesh,k_mesh,t),3);
% [lk_mesh,kl_mesh] = meshgrid(lam,k);
% loss_map(loss_map>1)=nan;
% figure;surf(lk_mesh,kl_mesh,loss_map);xlabel('$\lambda$');ylabel('k');zlabel('loss');xlim([min(lam),max(lam)]);ylim([min(k),max(k)]);
% hold on; plot(
%% weibull plot
% x = random(makedist('weibull',2,1.5),[1e4,1]);
% figure; histogram(x);
% [lam_mesh,k_mesh,x_mesh] = meshgrid(lam,k,x);
% loss_map = mean(loss(x_mesh,lam_mesh,k_mesh),3);
% loss_map(loss_map>1.5)=nan;
% figure;surf(lk_mesh,kl_mesh,loss_map);xlabel('$\lambda$');ylabel('k');zlabel('loss');xlim([min(lam),max(lam)]);ylim([min(k),max(k)]);
% x_trunc = x(x>1);
% [lam_mesh,k_mesh,x_mesh] = meshgrid(lam,k,x_trunc);
% loss_map = mean(loss_hat(x_mesh,lam_mesh,k_mesh,1),3);
% % loss_map(loss_map>1.5)=nan;
% figure;surf(lk_mesh,kl_mesh,loss_map);xlabel('$\lambda$');ylabel('k');zlabel('loss');xlim([min(lam),max(lam)]);ylim([min(k),max(k)]);
% % [xk_mesh,kx_mesh] = meshgrid(x,k);
% figure;opt_lambda = (sum(xk_mesh.^(kx_mesh-1),2)).^(1./k');plot(k,opt_lambda);
% figure;contour(loss)
%% weibull plot 2
m_ref = 1.5;
k_ref =2.4;
total_samples=1e2;
per_samples = .5;

f_ref = makedist('weibull',m_ref,k_ref);
x = random(f_ref,[total_samples,1]);
tr = icdf(f_ref,per_samples);
x_trunc= x(x>tr);
per_samples = 1- (length(x_trunc)./length(x));
m =.1:.05:10;
k = .1:.05:10;

[x_mesh,k_mesh,m_mesh]=meshgrid(x_trunc,k,m);
% loss_map = k_mesh*0;
% loss_hat_map = loss_map;
% loss_hat_reg_map = loss_map;
% loss_map_wrong = loss_map;
% reg_term = loss_map;

f = @(x,m,k) ((k./m).*((x./m).^(k-1)).*exp(-((x./m).^k)));
f_cdf = @(x,m,k) 1-exp(-(x./m).^k);
icdf_f = @(p,m,k) m.*(-log(1-p)).^(1./k);
p_hat = @(x,m,k,tr)(f(x,m,k)./(1-f_cdf(tr,m,k))).*((x>=tr)+eps);
reg = @(x,m,k,tr,p)(icdf_f(p,m,k)-tr).^2;
loss = @(x,m,k,tr,p,lambda)-squeeze(mean(log(p_hat(x,m,k,tr)),2))+lambda.*squeeze(mean(reg(x,m,k,tr,p),2));
loss_hat_reg_map = loss(x_mesh,m_mesh,k_mesh,tr,per_samples,.1);
loss_hat_map = loss(x_mesh,m_mesh,k_mesh,tr,per_samples,0);
loss_fun =  @(m,k)loss(x_trunc',m,k,tr,per_samples,0);
loss_reg_fun = @(m,k)loss(x_trunc',m,k,tr,per_samples,.1);
% for i_k = 1:length(k)
%     for i_m = 1:length(m)
%         f = makedist('weibull',m(i_m),k(i_k));
%         p = f.pdf(x);
%         p_hat = (f.pdf(x_trunc)./(1-f.cdf(tr))).*((x_trunc>=tr)+eps);
%         loss_map(i_k,i_m) = -mean(log(p));
%         loss_map_wrong(i_k,i_m) = -mean(log(f.pdf(x_trunc)));
%         loss_hat_map(i_k,i_m) = -mean(log(p_hat));
%         reg_term(i_k,i_m) = ((f.icdf(per_samples)-tr).^2);
%         loss_hat_reg_map(i_k,i_m) = -mean(log(p_hat))+5*reg_term(i_k,i_m);
%     end
% end
fun = @(theta)loss_fun(theta(1),theta(2));
options = optimoptions(@fmincon,'Display','off','Algorithm','interior-point','MaxFunctionEvaluations',10000);
[params,loss_val,~,out]=fmincon(fun,[1,1],[],[],[],[],[eps,eps],[],[],options);
fun = @(theta)loss_reg_fun(theta(1),theta(2));
options = optimoptions(@fmincon,'Display','off','Algorithm','interior-point','MaxFunctionEvaluations',10000);
[params_reg,loss_val,~,out]=fmincon(fun,[1,1],[],[],[],[],[eps,eps],[],[],options);
prctl_val = 10;
% loss_map(loss_map>prctile(loss_map(:),prctl_val))=nan;
loss_hat_map(loss_hat_map>prctile(loss_hat_map(:),prctl_val))=nan;
loss_hat_reg_map(loss_hat_reg_map>prctile(loss_hat_reg_map(:),prctl_val))=nan;
% loss_map_wrong(loss_map_wrong>prctile(loss_map_wrong(:),prctl_val))=nan;
% reg_term(reg_term>prctile(reg_term(:),prctl_val))=nan;
[m_mesh,k_mesh] = meshgrid(m,k);
view(0,90)
figure;
% subplot(2,2,1);surf(m_mesh,k_mesh,loss_map);title('Normal Loss');xlabel('$\lambda$');ylabel('k');zlabel('loss');xlim([min(m),max(m)]);ylim([min(k),max(k)]);
subplot(2,1,1);surf(m_mesh,k_mesh,loss_hat_map);view(0,90);xlabel('$\lambda$');title('Loss Hat');ylabel('k');zlabel('loss');xlim([min(m),max(m)]);ylim([min(k),max(k)]);hold on;
subplot(2,1,2);surf(m_mesh,k_mesh,loss_hat_reg_map);view(0,90);xlabel('$\lambda$');title('Loss with Regularization');ylabel('k');zlabel('loss');xlim([min(m),max(m)]);ylim([min(k),max(k)]);
% subplot(2,2,4);histogram(x);hold on;histogram(x_trunc);




% function loss = loss_fun(x,m,k,tr,per,lambda)
%     f = makedist('weibull',m(i_m),k(i_k));
%     p_hat = (f.pdf(x_trunc)./(1-f.cdf(tr))).*((x_trunc>=tr)+eps);
%     loss_hat_map = -mean(log(p_hat));
%     reg_term = ((f.icdf(per)-tr).^2);
%     loss_hat_reg_map= -mean(log(p_hat))+lambda*reg_term;
%     loss = loss_hat_reg;
% end
% subplot(2,3,5);surf(m_mesh,k_mesh,loss_map_wrong);xlabel('$\lambda$');title('Wrong Loss');ylabel('k');zlabel('loss');xlim([min(m),max(m)]);ylim([min(k),max(k)]);
% subplot(2,3,6);surf(m_mesh,k_mesh,reg_term);xlabel('$\lambda$');ylabel('k');title('Regularization');zlabel('loss');xlim([min(m),max(m)]);ylim([min(k),max(k)]);
