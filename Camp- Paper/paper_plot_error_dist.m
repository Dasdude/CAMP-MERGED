clc;close all;clear;
f = @(x,m,k) ((k./m).*((x./m).^(k-1)).*exp(-((x./m).^k)));
f_cdf = @(x,m,k) 1-exp(-(x./m).^k);
icdf_f = @(p,m,k) m.*(-log(1-p)).^(1./k);
p_hat = @(x,m,k,tr)(f(x,m,k)./(1-f_cdf(tr,m,k))).*((x>=tr)+eps);
reg = @(x,m,k,tr,p)(icdf_f(p,m,k)-tr).^2;
loss = @(x,m,k,tr,p,lambda)-squeeze(mean(log(p_hat(x,m,k,tr)),2))+lambda.*squeeze(mean(reg(x,m,k,tr,p),2));
% loss_fun =  @(x,m,k)loss(x_trunc',m,k,tr,per_samples,0);
% loss_reg_fun = @(x,m,k)loss(x_trunc',m,k,tr,per_samples,.1);
ent = @(m,k)eulergamma*(1-(1./k)) - log(m./k)+1;


m_ref = 1.5;
k_ref =2.4;
f_ref = makedist('weibull',m_ref,k_ref);
total_samples=1e2;
per_samples_list = 0:.3:1;
total_trials = 10;
ent_ref = ent(m_ref,k_ref);
err_reg  = zeros(length(per_samples_list),total_trials)*nan;
err_hat = err_reg;
err_new = err_reg;
for per_samples_idx = 1:length(per_samples_list)
    per_samples = per_samples_list(per_samples_idx);
    sprintf('PER: %d',per_samples)
    for i = 1:total_trials
        x = random(f_ref,[total_samples,1]);
        tr = icdf(f_ref,per_samples);
        x_trunc= x(x>tr);
        if length(x_trunc)<5
            continue
        end
        per = 1-(length(x_trunc)./length(x));
        if per==1
            continue
        end
        loss_fun =  @(x,m,k)loss(x,m,k,tr,per,0);
        loss_reg_fun = @(x,m,k)loss(x,m,k,tr,per,.1);
        loss_new_fun = @(x,m,k)categorical_loss(x,per_samples,m,k);
        fun = @(theta)loss_fun(x_trunc',theta(1),theta(2));
        options = optimoptions(@fmincon,'Display','off','Algorithm','interior-point','MaxFunctionEvaluations',10000);
        [params,loss_val,~,out]=fmincon(fun,[1,1],[],[],[],[],[eps,eps],[],[],options);
        fun = @(theta)loss_reg_fun(x_trunc',theta(1),theta(2));
        options = optimoptions(@fmincon,'Display','off','Algorithm','interior-point','MaxFunctionEvaluations',10000);
        [params_reg,loss_val_reg,~,out]=fmincon(fun,[1,1],[],[],[],[],[eps,eps],[],[],options);
        
        fun = @(theta)loss_new_fun(x_trunc',theta(1),theta(2));
        options = optimoptions(@fmincon,'Display','off','Algorithm','interior-point','MaxFunctionEvaluations',10000);
        [params_new,loss_val_new,~,out]=fmincon(fun,[1,1],[],[],[],[],[eps,eps],[],[],options);
        loss_val_no_reg = loss_val;
        loss_val_reg = loss_fun(x_trunc',params_reg(1),params_reg(2));
        ll_ref = loss_fun(x_trunc',m_ref,k_ref);
        
        err_reg(per_samples_idx,i) = sqrt(sum((params_reg-[m_ref,k_ref]).^2));
        err_hat(per_samples_idx,i) = sqrt(sum((params-[m_ref,k_ref]).^2));
        err_new(per_samples_idx,i) = sqrt(sum((params_new-[m_ref,k_ref]).^2));
       
    end
     
     sprintf('mean reg: %d hat: %d',mean(err_reg(per_samples_idx,:),'omitnan'),mean(err_hat(per_samples_idx,:),'omitnan'))
%      sprintf('std reg: %d hat: %d',std(err_hat(per_samples_idx,:),'omitnan'),std(err_hat(per_samples_idx,:),'omitnan'))
end
err_reg_mean = mean(err_reg,2,'omitnan');
err_hat_mean = mean(err_hat,2,'omitnan');
err_new_mean = mean(err_new,2,'omitnan');
figure;plot(per_samples_list,[err_reg_mean,err_hat_mean,err_new_mean]);legend('reg','hat','new');
err_reg_std = std(err_reg','omitnan');
err_hat_std = std(err_hat','omitnan');
tmp = per_samples_list'.*ones(1,total_trials);
tmp = tmp(:);
figure;subplot(3,1,1);scatter(tmp,err_reg(:),20);subplot(3,1,2);scatter(tmp,err_hat(:),20);subplot(3,1,3);scatter(tmp,err_new(:),20);
% figure;plot(per_samples_list,[err_hat_std',err_reg_std']);