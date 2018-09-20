clc
close all
clear
addpath(genpath('.'))
%% File Names
PLOT_EDITED_PARAMS = 0;
experiment_name = 'Aug 16 General';
addpath(genpath('.'))
%% Distribution
dist_obj_nak = distribution_type_class(@(x)makedist('nakagami',x(1),x(2)),'nakagami',{'mu','omega'},[0.5,0],[inf,inf]);
dist_obj_logn = distribution_type_class(@(x)makedist('lognormal',x(1),x(2)),'lognormal',{'mu','sigma'},[-inf,0],[3,3]);
dist_obj_wei = distribution_type_class(@(x)makedist('weibull',x(1),x(2)),'weibull',{'A','B'},[0,0],[5,5]);
dist_obj_ri = distribution_type_class(@(x)makedist('rician',x(1),x(2)),'rician',{'s','sigma'},[0,0],[5,5]);
dist_obj_ray = distribution_type_class(@(x)makedist('rayleigh',x(1)),'rayleigh',{'B'},[0],[5]);
dist_obj_cell = {dist_obj_nak,dist_obj_logn,dist_obj_wei,dist_obj_ri,dist_obj_ray};
dist_obj = dist_obj_cell{5};
%% Constant Iniit
same_low_up = 15;
same_med_up = 30;

opposite_low_up = 15;
opposite_med_up = 30;

SAME_DENS_LOW = {'same','low','0',same_low_up,1.0043,2.0108};
SAME_DENS_MED = {'same','med',same_low_up,same_med_up,1.0091,2.0237};
SAME_DENS_HIGH = {'same','high',same_med_up,100,1.0036,2.0576};
OP_DENS_LOW = {'opposite','low','0',opposite_low_up,1,2.1434};
OP_DENS_MED = {'opposite','med',opposite_low_up,opposite_med_up,1,2.1904};
OP_DENS_HIGH = {'opposite','high',opposite_med_up,100,1,2.26};
mode_list = {SAME_DENS_LOW,SAME_DENS_MED,SAME_DENS_HIGH,OP_DENS_LOW,OP_DENS_MED,OP_DENS_HIGH};
%% Distribution
% dist_obj_nak = distribution_type_class(@(x)makedist('nakagami',x(1),x(2)),'nakagami',{'mu','omega'},[0.5,0],[inf,inf]);
% dist_obj_logn = distribution_type_class(@(x)makedist('lognormal',x(1),x(2)),'lognormal',{'mu','sigma'},[-inf,0],[3,3]);
% dist_obj_wei = distribution_type_class(@(x)makedist('weibull',x(1),x(2)),'weibull',{'A','B'},[0,0],[5,5]);
% dist_obj_ri = distribution_type_class(@(x)makedist('rician',x(1),x(2)),'rician',{'s','sigma'},[0,0],[5,5]);
% dist_obj_ray = distribution_type_class(@(x)makedist('rayleigh',x(1)),'rayleigh',{'B'},[0],[5]);
% dist_obj_cell = {dist_obj_nak,dist_obj_logn,dist_obj_wei,dist_obj_ri,dist_obj_ray};
if PLOT_EDITED_PARAMS
    output_extension = '_edit.png';
    param_file_name = 'Parameters_edit.mat';
else
    output_extension = '.png';
    param_file_name = 'Parameters.mat';
end
noise_level = -98;
pkt_size = -1;
TRUNCATION_VALUE = -90;
censor_function_handle = @(x)censor_function(x,noise_level,pkt_size,TRUNCATION_VALUE);
for mode_index = 1:6
    mode = mode_list{mode_index};
    mode_name = sprintf('%s',mode{1});
    %% Dataset prepare
        display('Data Prepare Phase')
        dataset_file_path  = sprintf('%s/s %d %d o %d %d/data results/%s %s.csv','Seperated DensityPER',same_low_up,same_med_up,opposite_low_up,opposite_med_up,mode{1},mode{2});;
        
        
        csv_data = readtable(dataset_file_path,'ReadVariableNames',true);
    for dist_index = 1:5       
        close all
        dist_obj = dist_obj_cell{dist_index};
        dist_name = dist_obj.dist_name;
        minimal_experiment_name = [replace(mode{1},'_',' '),' ',dist_obj.dist_name];
        relative_experiment_folder_path = fullfile('Plots',experiment_name,mode{1},mode{2},dist_obj.dist_name);
        parameter_folder = fullfile(relative_experiment_folder_path,'Results');
        
        parameter_path = fullfile(parameter_folder,param_file_name);
        plot_folder_path = relative_experiment_folder_path;
        %% Load Params
        load(parameter_path);
        %% Parameters
        d_max = length(fading_params);
        TRUNCATION_VALUE=-90;
        LIGHT_SPEED=3*10^8;
        lambda=LIGHT_SPEED/CARRIER_FREQ;
        dataset_mat_dirty = [csv_data.Range,csv_data.RSS];
        dataset_mat_dirty(dataset_mat_dirty(:,2)==-101,2) = -999;
        dataset_cell_dirty = data_mat_cell(dataset_mat_dirty,d_max);
        [dataset_cell,per,packet_loss_stat]=censor_data(dataset_cell_dirty,censor_function_handle);
        data_dbm_cell = dataset_cell;
        data_dbm_cell = data_dbm_cell(1:d_max);
        data_dbm_mean = funoncellarray1input(data_dbm_cell,@mean);
        data_dbm_std = funoncellarray1input(data_dbm_cell,@std);

        %% Pathloss
        pathloss = pathloss_gen_2ray(TX_HEIGHT,RX_HEIGHT,EPSILON,ALPHA,lambda,d_max);
%         pathloss = 1.3*(20*log10(1:d_max)+20*log10(CARRIER_FREQ)+20*log10(4*pi/LIGHT_SPEED));
        %% Extract Fading
        data_fading_dbm = extract_fading(data_dbm_cell,pathloss,TX_POWER);
        %% Generate Data
        generated_fading_linear = sample_generator(dist_obj,fading_params,1e3);
        generated_fading_dbm = linear2dbm(generated_fading_linear);
        generated_rssi_dbm = add_fading(pathloss,generated_fading_dbm,TX_POWER);
        [generated_rssi_dbm_truncated,generated_per,gen_pl_stat] = censor_data(generated_rssi_dbm,censor_function_handle);
        generated_rssi_dbm_mean = funoncellarray1input(generated_rssi_dbm,@mean);
        
           %% Pathloss Compare Plot
        figure;subplot(2,1,1);plot(generated_rssi_dbm_mean);hold;plot(data_dbm_mean);title([minimal_experiment_name,'Mean Comparison']);grid on;legend('Model','Field');subplot(2,1,2);plot(aprx_per);title('PER');ylabel('RSS');saveas(gcf,[plot_folder_path,'/','Mean Model vs Field',output_extension]);
        %% Percentile Plot
        non_trunc_ylim = [-130,-30];
        percentiles_generated = percentile_array([5,10,25,50,75,90,95],generated_rssi_dbm);
        percentiles_generated_trunc = percentile_array([5,10,25,50,75,90,95],generated_rssi_dbm_truncated);
        percentiles_rssi = percentile_array([5,10,25,50,75,90,95],data_dbm_cell);
        percentiles_rssi_per_inc = percentile_array_per([5,10,25,50,75,90,95],data_dbm_cell,per*100);
%         percentiles_rssi_gen_per_inc = percentile_array_per([5,10,25,50,75,90,95],generated_rssi_dbm_truncated,generated_per*100);
        figure;plot(percentiles_generated(:,[2,4,6]));hold on ;plot(percentiles_rssi_per_inc(:,[2,4,6]));grid on;ylim(non_trunc_ylim);title([minimal_experiment_name,'Percentile']);ylabel('RSS (dbm)');xlabel('Distance (m)');legend('10% model','50% model','90% model','10% field','50% field','90% field');saveas(gcf,[plot_folder_path,'/','Percentile RSSI 10',output_extension]);
        figure;plot(percentiles_generated_trunc(:,[2,4,6]));hold on ;plot(percentiles_rssi(:,[2,4,6]));grid on;title([minimal_experiment_name,'Truncated Percentile']);ylabel('RSS (dbm)');xlabel('Distance (m)');legend('10% model','50% model','90% model','10% field','50% field','90% field');saveas(gcf,[plot_folder_path,'/','Percentile RSSI Truncated 10',output_extension]);
        figure;plot(percentiles_generated(:,[3,4,5]));hold on ;plot(percentiles_rssi_per_inc(:,[3,4,5]));grid on;ylim(non_trunc_ylim);title([minimal_experiment_name,'Percentile']);ylabel('RSS (dbm)');xlabel('Distance (m)');legend('25% model','50% model','75% model','25% field','50% field','75% field');saveas(gcf,[plot_folder_path,'/','Percentile RSSI 25',output_extension]);
        figure;plot(percentiles_generated_trunc(:,[3,4,5]));hold on ;plot(percentiles_rssi(:,[3,4,5]));grid on;title([minimal_experiment_name,'Truncated Percentile']);ylabel('RSS (dbm)');xlabel('Distance (m)');legend('25% model','50% model','75% model','25% field','50% field','75% field');saveas(gcf,[plot_folder_path,'/','Percentile RSSI Truncated 25',output_extension]);
        figure;plot(percentiles_generated(:,[1,4,7]));hold on ;plot(percentiles_rssi_per_inc(:,[1,4,7]));grid on;ylim(non_trunc_ylim);title([minimal_experiment_name,'Percentile']);ylabel('RSS (dbm)');xlabel('Distance (m)');legend('5% model','50% model','95% model','5% field','50% field','95% field');saveas(gcf,[plot_folder_path,'/','Percentile RSSI 5',output_extension]);
        figure;plot(percentiles_generated_trunc(:,[1,4,7]));hold on ;plot(percentiles_rssi(:,[1,4,7]));grid on;title([minimal_experiment_name,'Truncated Percentile']);ylabel('RSS (dbm)');xlabel('Distance (m)');legend('5% model','50% model','95% model','5% field','50% field','95% field');saveas(gcf,[plot_folder_path,'/','Percentile RSSI Truncated 5',output_extension]);
        %% PER Plot
        figure;plot(packet_loss_stat(:,2));hold;plot(packet_loss_stat(:,2)-packet_loss_stat(:,1));xlabel('Distance(m)');ylabel('Number of Samples');grid on;title([minimal_experiment_name,'Total Samples vs Received Samples']);legend('Total Samples','Recieved Samples');saveas(gcf,[plot_folder_path,'/','Samples Received vs Total',output_extension]);
        figure; plot(100*generated_per);hold on;plot(100*per);grid on;title([minimal_experiment_name,'PER Value Comparison']);ylabel('Rate');xlabel('Distance (m)');legend('Model','Field','Location','northwest');saveas(gcf,[plot_folder_path,'/','PER Comparison',output_extension]);
        figure;plot(loss_vals);title([minimal_experiment_name,'loss']);grid on;saveas(gcf,[plot_folder_path,'/','Loss',output_extension]);
        %% Plot Parameters
        for param_idx = 1:dist_obj.get_dof
            param_name = dist_obj.dist_params_names{param_idx};
            figure;plot(fading_params(:,param_idx));title(sprintf('%s Parameter: %s - Distance',minimal_experiment_name,param_name));grid on;xlabel('Distance (m)');ylabel(sprintf('%s Value',param_name));saveas(gcf,fullfile(plot_folder_path,sprintf('%s_distance%s',param_name,output_extension)));
        end
        
        %% Plot Normalized Likelihood
%         loglikelihood_set([1,1],data_dbm_cell(100),.9,dist_obj.dist_handle)
%         ll_fun_handle = @(x)loglikelihood_set(,data_dbm_cell(100),.9,dist_obj.dist_handle)
%         fading_trunc_val_dbm = -94+pathloss-TX_POWER;
%         tmp_input_cell_array = {data_fading_dbm,fading_params,-inf*ones(length(generated_fading_dbm),1)};
%         ll_fading_dbm = funonarray(ll_fun_handle,tmp_input_cell_array);
%         tmp_input_cell_array = {data_fading_dbm,fading_params,fading_trunc_val_dbm'};
%         ll_fading_dbm_truncated = funonarray(ll_fun_handle,tmp_input_cell_array);
%         % Plot
%         figure;plot(ll_fading_dbm);hold on ;plot(ll_fading_dbm_truncated);grid on;title([minimal_experiment_name,'Log Likelihood RSSI - Distance']);legend('Full Distribution','Truncated Distribution');saveas(gcf,[plot_folder_path,'/','LL RSSI',output_extension]);
%         %% KS-Test
%         h = funoncellarray2input(generated_rssi_dbm_truncated,data_dbm_cell,@(x,y)kstest2(x,y,'Alpha',.01));
%         figure;plot(h);title('KS-Test');xlabel('Distance (m)');ylim([-1,2]);ylabel('Null Hyphotesis State');grid on;saveas(gcf,[plot_folder_path,'/','KS-Test',output_extension]);
    end
end