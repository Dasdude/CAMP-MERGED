clc
close all
clear
addpath(genpath('.'))
%% File Names
% mode_index =6;
pkt_size = -1;
noise_level = -98;
for mode_index = 1:6
    close all
same_low_up = 15;
same_med_up = 30;

opposite_low_up = 15;
opposite_med_up = 30;

SAME_DENS_LOW = {'same','low',0,same_low_up,1.0043,2.0108};
SAME_DENS_MED = {'same','med',same_low_up,same_med_up,1.0091,2.0237};
SAME_DENS_HIGH = {'same','high',same_med_up,100,1.0036,2.0576};
OP_DENS_LOW = {'opposite','low',0,opposite_low_up,1,2.1434};
OP_DENS_MED = {'opposite','med',opposite_low_up,opposite_med_up,1,2.1904};
OP_DENS_HIGH = {'opposite','high',opposite_med_up,100,1,2.26};
mode_list = {SAME_DENS_LOW,SAME_DENS_MED,SAME_DENS_HIGH,OP_DENS_LOW,OP_DENS_MED,OP_DENS_HIGH};
mode = mode_list{mode_index};
experiment_name = 'Camp Highway Model Final July 30- 14dbm';
minimal_experiment_name = [upper(mode{1}),' ',upper(mode{2})];
mode_name = [mode{1},' Direction ',mode{2},' Density ',num2str(mode{3}),' to ',num2str(mode{4})];
parameter_folder = ['Plots/',experiment_name,'/',mode_name,'/Results'];
parameter_path = [parameter_folder,'/Parameters.mat'];
%% Load Params
load(parameter_path);
%% Parameters
d_max = 800;
TRUNCATION_VALUE=-93;
LIGHT_SPEED=3*10^8;
lambda=LIGHT_SPEED/CARRIER_FREQ;
%% Dataset prepare
display(['Data Prepare Phase ',mode_name]);
file_string = sprintf('%s/s %d %d o %d %d/data results/%s %s.csv','Seperated DensityPER',same_low_up,same_med_up,opposite_low_up,opposite_med_up,mode{1},mode{2});
input  = file_string;
file_name_string = sprintf('%s/%s Direction %s Density %d to %d',experiment_name,mode{1},mode{2},mode{3},mode{4});
csv_data = readtable(input,'ReadVariableNames',true);
dataset_mat_dirty = [csv_data.Range,csv_data.RSS];
any(isnan(dataset_mat_dirty))

any(dataset_mat_dirty(:)<-100)
dataset_cell_dirty = data_mat_cell(dataset_mat_dirty,d_max);
packet_loss_stat = per_calc(dataset_cell_dirty,-100);
% per = packet_loss_stat(:,1)./packet_loss_stat(:,2)*100;
censor_function_handle = @(x)censor_function(x,noise_level,pkt_size,TRUNCATION_VALUE);
[dataset_cell,per] = censor_data(dataset_cell_dirty,censor_function_handle);
per = per*100;
data_dbm_cell = dataset_cell;
data_dbm_cell = data_dbm_cell(1:d_max);
data_dbm_mean = funoncellarray1input(data_dbm_cell,@mean);
data_dbm_std = funoncellarray1input(data_dbm_cell,@std);
%% Check Params
mu_correct = all(nakagami_mu == fading_params(:,1));
omega_correct = all(nakagami_omega == fading_params(:,2));
if ~(mu_correct&omega_correct)
    error('Parameters Dont Match')
end
%% Pathloss
pathloss = pathloss_gen_2ray(TX_HEIGHT,RX_HEIGHT,EPSILON,ALPHA,lambda,d_max);
%% Extract Fading
TX_POWER=14;
data_fading_dbm = extract_fading(data_dbm_cell,pathloss,TX_POWER);
%% Generate Data
fading_params(:,2) = fading_params(:,2)*10^0.3; % the parameters stored are for 17 dbm therefore for 14 dbm the parameters should be changed
generated_fading_linear = nakagami_generator(fading_params,1e3);
generated_fading_dbm = linear2dbm(generated_fading_linear);
generated_rssi_dbm = add_fading(pathloss,generated_fading_dbm,TX_POWER);
% generated_rssi_dbm_truncated = truncate_data_cell(generated_rssi_dbm,TRUNCATION_VALUE);
[generated_rssi_dbm_truncated,generated_per] = censor_data(generated_rssi_dbm,censor_function_handle);
generated_per = generated_per*100;
generated_rssi_dbm_mean = funoncellarray1input(generated_rssi_dbm,@mean);
generated_total_samples = funoncellarray1input(generated_rssi_dbm,@length);
generated_received_samples = funoncellarray1input(generated_rssi_dbm_truncated,@length);
% generated_recieved_samples = 
% generated_per = (1-(generated_received_samples./generated_total_samples))*100;
%% Pathloss Compare Plot
figure;subplot(2,1,1);plot(generated_rssi_dbm_mean);hold;plot(data_dbm_mean);title([ minimal_experiment_name,' Mean Comparison']);grid on;legend('Model','Field');subplot(2,1,2);plot(generated_per);hold;plot(per);ylim([0,100]);legend('Model','Field','Location','southeast');title([ minimal_experiment_name,' PER Comparison']);xlabel('TxRx Distance(m)');ylabel('PER(%)');grid on;saveas(gcf,['Plots/',file_name_string,'/','Mean Model vs Field.png']);
figure;plot(data_dbm_mean);hold;plot(TX_POWER-pathloss);xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');legend('Field Data Mean','Estimated Pathloss');title([ minimal_experiment_name,' Estimated Pathloss vs Mean RSS of Field Samples']);grid on;saveas(gcf,['Plots/',file_name_string,'/','Pathloss vs Field.png']);
%% Percentile Plot
percentiles_generated = percentile_array([5,10,25,50,75,90,95],generated_rssi_dbm);
percentiles_generated_trunc = percentile_array([5,10,25,50,75,90,95],generated_rssi_dbm_truncated);
percentiles_rssi = percentile_array([5,10,25,50,75,90,95],data_dbm_cell);
figure;plot(percentiles_generated(:,[1,4,7]));hold on ;plot(percentiles_rssi(:,[1,4,7]));grid on;title([ minimal_experiment_name,' RSS Percentile']);xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');legend('5% model','50% model','95% model','5% field','50% field','95% field');saveas(gcf,['Plots/',file_name_string,'/','Percentile RSS 5.png']);
figure;plot(percentiles_generated_trunc(:,[1,4,7]));hold on ;plot(percentiles_rssi(:,[1,4,7]));grid on;title([ minimal_experiment_name,' Truncated RSS Percentile']);xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');legend('5% model','50% model','95% model','5% field','50% field','95% field');saveas(gcf,['Plots/',file_name_string,'/','Percentile RSS Truncated 5.png']);
figure;plot(percentiles_generated(:,[2,4,6]));hold on ;plot(percentiles_rssi(:,[2,4,6]));grid on;title([ minimal_experiment_name,' RSS Percentile']);xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');legend('10% model','50% model','90% model','10% field','50% field','90% field');saveas(gcf,['Plots/',file_name_string,'/','Percentile RSS 10.png']);
figure;plot(percentiles_generated_trunc(:,[2,4,6]));hold on ;plot(percentiles_rssi(:,[2,4,6]));grid on;title([ minimal_experiment_name,' Truncated RSS Percentile']);xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');legend('10% model','50% model','90% model','10% field','50% field','90% field');saveas(gcf,['Plots/',file_name_string,'/','Percentile RSS Truncated 10.png']);
figure;plot(percentiles_generated(:,[3,4,5]));hold on ;plot(percentiles_rssi(:,[3,4,5]));grid on;title([ minimal_experiment_name,'  RSS Percentile']);xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');legend('25% model','50% model','75% model','25% field','50% field','75% field');saveas(gcf,['Plots/',file_name_string,'/','Percentile RSS 25.png']);
figure;plot(percentiles_generated_trunc(:,[3,4,5]));hold on ;plot(percentiles_rssi(:,[3,4,5]));grid on;title([ minimal_experiment_name,' Truncated RSS Percentile']);xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');legend('25% model','50% model','75% model','25% field','50% field','75% field');saveas(gcf,['Plots/',file_name_string,'/','Percentile RSS Truncated 25.png']);
%% PER Plot
figure;plot(packet_loss_stat(:,2));hold;plot(packet_loss_stat(:,2)-packet_loss_stat(:,1));grid on;title([ minimal_experiment_name,' Total Samples vs Received Samples']);legend('Total Samples','Recieved Samples');saveas(gcf,['Plots/',file_name_string,'/','Samples Received vs Total.png']);
figure; plot(generated_per);hold on;plot(per);grid on;title([ minimal_experiment_name,' PER Value Comparison']);ylim([0,100]);ylabel('PER(%)');xlabel('Range');legend('Model','Field','Location','northwest');saveas(gcf,['Plots/',file_name_string,'/','PER Comparison.png']);

%% Scatter Plot
generated_rssi_dbm_mat = data_cell2mat(generated_rssi_dbm);
generated_rssi_dbm_truncated_mat = data_cell2mat(generated_rssi_dbm_truncated);
field_rssi_dbm_mat = data_cell2mat(data_dbm_cell);
figure;scatter(generated_rssi_dbm_mat(:,1),generated_rssi_dbm_mat(:,2),1.5,'.');grid on;title([minimal_experiment_name, 'Generated RSS From Model']);xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');saveas(gcf,['Plots/',file_name_string,'/','scatter generated.png']);
figure;scatter(generated_rssi_dbm_truncated_mat(:,1),generated_rssi_dbm_truncated_mat(:,2),1.5,'.');grid on;title([ minimal_experiment_name,' Generated Censored RSS From Model']);xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');ylim([min(generated_rssi_dbm_truncated_mat(:,2)),-20]);saveas(gcf,['Plots/',file_name_string,'/','scatter generated trunc.png']);
figure;scatter(field_rssi_dbm_mat(:,1),field_rssi_dbm_mat(:,2),1.5,'.');grid on;xlabel('Range(m)');ylabel('RSS(dbm)');title([ minimal_experiment_name,' Field Data']);xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');ylim([min(field_rssi_dbm_mat(:,2)),-20]);saveas(gcf,['Plots/',file_name_string,'/','scatter field Discrete.png']);
figure;scatter(field_rssi_dbm_mat(:,1),field_rssi_dbm_mat(:,2)+rand(size(field_rssi_dbm_mat(:,2)))-.5,1.5,'.');grid on;xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');title([ minimal_experiment_name,' Field RSS + Added Uniform Noise~[-.5,.5])']);ylim([min(field_rssi_dbm_mat(:,2)),-20]);saveas(gcf,['Plots/',file_name_string,'/','scatter field uniform noise.png']);
%% Plot Nakagami Parameters
figure;plot(fading_params(:,1));title([ minimal_experiment_name,' Mu - Distance']);xlabel('TxRx Distance(m)');ylabel('Mu Value');grid on;saveas(gcf,['Plots/',file_name_string,'/','mu_distance.png']);
figure;plot(fading_params(:,2));title([ minimal_experiment_name,' Omega - Distance']);grid on;xlabel('TxRx Distance(m)');ylabel('Omega Value');saveas(gcf,['Plots/',file_name_string,'/','Omega_distance.png']);
%% Plot Normalized Likelihood
figure;plot(loss_vals);title([ minimal_experiment_name,' Loss Function']);xlabel('TxRx Distance(m)');grid on;saveas(gcf,['Plots/',file_name_string,'/','Loss.png']);
ll_fun_handle = @(x)loglikelihood_samples(x{1},'lognakagami',x{2},x{3});
fading_trunc_val_dbm = -94+pathloss-TX_POWER;
tmp_input_cell_array = {data_fading_dbm,fading_params,-inf*ones(length(generated_fading_dbm),1)};
ll_fading_dbm = funonarray(ll_fun_handle,tmp_input_cell_array);
tmp_input_cell_array = {data_fading_dbm,fading_params,fading_trunc_val_dbm'};
ll_fading_dbm_truncated = funonarray(ll_fun_handle,tmp_input_cell_array);
% Plot
figure;plot(ll_fading_dbm);hold on ;plot(ll_fading_dbm_truncated);grid on;xlabel('TxRx Distance(m)');title([ minimal_experiment_name,' Log Likelihood']);legend('Full Distribution','Truncated Distribution');saveas(gcf,['Plots/',file_name_string,'/','LL RSS.png']);

save(['Plots/',file_name_string,'/','data.mat'],'percentiles_generated','percentiles_generated_trunc','per','generated_per','data_dbm_cell','percentiles_rssi','generated_rssi_dbm');

%% Parameters 
nakagami_mu = fading_params(:,1);
nakagami_omega = fading_params(:,2);

save(['Plots/',file_name_string,'/Results/','Parameters_14dbm.mat'],'TX_HEIGHT','RX_HEIGHT','tworay_pathloss_alpha','tworay_pathloss_epsilon','TX_POWER','CARRIER_FREQ','nakagami_mu','nakagami_omega')
end