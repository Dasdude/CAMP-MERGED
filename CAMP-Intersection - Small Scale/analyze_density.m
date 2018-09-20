clc
close all
clear
addpath(genpath('.'))
%% File Names
experiment_name = 'AUG13';
SAME_LEG = {'same_leg'};
DIF_LEG_LOS = {'dif_leg_los'};
DIF_LEG_NLOS = {'dif_leg_nlos'};
mode_list = {SAME_LEG,DIF_LEG_LOS,DIF_LEG_NLOS};

%% Distribution
dist_obj_nak = distribution_type_class(@(x)makedist('nakagami',x(1),x(2)),'nakagami',{'mu','omega'},[0.5,0],[inf,inf]);
dist_obj_logn = distribution_type_class(@(x)makedist('lognormal',x(1),x(2)),'lognormal',{'mu','sigma'},[-inf,0],[3,3]);
dist_obj_wei = distribution_type_class(@(x)makedist('weibull',x(1),x(2)),'weibull',{'A','B'},[0,0],[5,5]);
dist_obj_ri = distribution_type_class(@(x)makedist('rician',x(1),x(2)),'rician',{'s','sigma'},[0,0],[5,5]);
dist_obj_ray = distribution_type_class(@(x)makedist('rayleigh',x(1)),'rayleigh',{'B'},[0],[5]);
dist_obj_cell = {dist_obj_nak,dist_obj_logn,dist_obj_wei,dist_obj_ri,dist_obj_ray};

noise_level = -98;
pkt_size = -1;
TRUNCATION_VALUE = -90;
censor_function_handle = @(x)censor_function(x,noise_level,pkt_size,TRUNCATION_VALUE);
time_quant_val = 2;
for mode_index = 1:3
    mode = mode_list{mode_index};
    mode_name = sprintf('%s',mode{1});
    minimal_experiment_name = [replace(mode{1},'_',' ')];
    plot_folder_path = fullfile('Plots','DataAnalysis','Density Analysis',mode{1});
    mkdir(plot_folder_path);
    %% Dataset prepare
    display('Data Prepare Phase')
    dataset_file_path = sprintf('Dataset/Ehsan/%s.csv',mode{1});
    csv_data = readtable(dataset_file_path,'ReadVariableNames',true);
    a=datevec(csv_data.Timestamp/60/60/24/1000) + [1970 0 0 0 0 0];
    csv_data.Hour = a(:,4);
    
    csv_data.twoHour = floor(a(:,4)/time_quant_val)*time_quant_val;
    time_set = unique(csv_data.twoHour);
    d_max = floor(prctile(csv_data.TxRxDistance(csv_data.RSS<300),99.99));
    time_name_cell_array = {};
    per = zeros(d_max,length(time_set));
    percentiles_rssi=zeros(d_max,7,length(time_set)) ;
    percentiles_rssi_per_inc=zeros(d_max,7,length(time_set)) ;
    for time_idx = 1:length(time_set)
        
        time_val = time_set(time_idx);
        time_name_cell_array{time_idx} = sprintf('Time %g:%g',time_val,mod(time_val+time_quant_val,24)) ;
        data_sel = csv_data(csv_data.twoHour==time_val,:);
        dataset_mat_dirty = [data_sel.TxRxDistance,data_sel.RSS];
        dataset_mat_dirty(dataset_mat_dirty(:,2)>300,2) = -999;
        dataset_mat_dirty(dataset_mat_dirty(:,2)<-130,2) = -999;
        dataset_cell_dirty = data_mat_cell(dataset_mat_dirty,d_max);
        [dataset_cell,per(:,time_idx),packet_loss_stat]=censor_data(dataset_cell_dirty,censor_function_handle);
        data_dbm_cell = dataset_cell;
        data_dbm_cell = data_dbm_cell(1:d_max);
        data_dbm_mean = funoncellarray1input(data_dbm_cell,@mean);
        data_dbm_std = funoncellarray1input(data_dbm_cell,@std);
%         d_max = min(d_max,length(fading_params));
        close all     
        %% Pathloss Compare Plot
%         figure;subplot(2,1,1);plot(generated_rssi_dbm_mean);hold;plot(data_dbm_mean);title([minimal_experiment_name,'Mean Comparison']);grid on;legend('Model','Field');subplot(2,1,2);plot(aprx_per);title('PER');ylabel('RSS');saveas(gcf,[plot_folder_path,'/','Mean Model vs Field.png']);
        %% Percentile Plot
        non_trunc_ylim = [-130,-30];
        percentiles_rssi(:,:,time_idx) = percentile_array([5,10,25,50,75,90,95],data_dbm_cell);
        percentiles_rssi_per_inc(:,:,time_idx) = percentile_array_per([5,10,25,50,75,90,95],data_dbm_cell,per*100);
    end
        figure; plot(squeeze(percentiles_rssi_per_inc(:,4,:)));grid on;title([minimal_experiment_name,'Median']);ylabel('RSS (dbm)');xlabel('Distance (m)');legend(time_name_cell_array);saveas(gcf,[plot_folder_path,'/','Median Comparison.png']);
        figure; plot(squeeze(percentiles_rssi_per_inc(:,6,:)));grid on;title([minimal_experiment_name,'Time Based 90th Percentile Including Lost Packets']);ylabel('RSS (dbm)');xlabel('Distance (m)');legend(time_name_cell_array);saveas(gcf,[plot_folder_path,'/','90th Comparison.png']);
        figure; plot(squeeze(percentiles_rssi_per_inc(:,2,:)));grid on;title([minimal_experiment_name,'Time Based 10th Percentile Including Lost Packets']);ylabel('RSS (dbm)');xlabel('Distance (m)');legend(time_name_cell_array);saveas(gcf,[plot_folder_path,'/','10th Comparison.png']);
        
        figure; plot(squeeze(percentiles_rssi(:,4,:)));grid on;title([minimal_experiment_name,'Median']);ylabel('RSS (dbm)');xlabel('Distance (m)');legend(time_name_cell_array);saveas(gcf,[plot_folder_path,'/','PLOSS Median Comparison.png']);
        figure; plot(squeeze(percentiles_rssi(:,6,:)));grid on;title([minimal_experiment_name,'Time Based 90th Percentile']);ylabel('RSS (dbm)');xlabel('Distance (m)');legend(time_name_cell_array);saveas(gcf,[plot_folder_path,'/','PLOSS 90th Comparison.png']);
        figure; plot(squeeze(percentiles_rssi(:,2,:)));grid on;title([minimal_experiment_name,'Time Based 10th Percentile']);ylabel('RSS (dbm)');xlabel('Distance (m)');legend(time_name_cell_array);saveas(gcf,[plot_folder_path,'/','PLOSS 10th Comparison.png']);
        figure; plot(per*100);grid on;title([minimal_experiment_name,'PER']);ylabel('PER (%)');xlabel('Distance (m)');legend(time_name_cell_array);saveas(gcf,[plot_folder_path,'/','PER.png']);
        figure; plot(medfilt1(per*100,10));grid on;title([minimal_experiment_name,' Smooth PER']);ylabel('PER (%)');xlabel('Distance (m)');legend(time_name_cell_array);saveas(gcf,[plot_folder_path,'/','PER Smooth.png']);
        
%         percentiles_rssi_gen_per_inc = percentile_array_per([5,10,25,50,75,90,95],generated_rssi_dbm_truncated,generated_per*100);
       
        %% Plot Normalized Likelihood
%         loglikelihood_set([1,1],data_dbm_cell(100),.9,dist_obj.dist_handle)
%         ll_fun_handle = @(x)loglikelihood_set(,data_dbm_cell(100),.9,dist_obj.dist_handle)
%         fading_trunc_val_dbm = -94+pathloss-TX_POWER;
%         tmp_input_cell_array = {data_fading_dbm,fading_params,-inf*ones(length(generated_fading_dbm),1)};
%         ll_fading_dbm = funonarray(ll_fun_handle,tmp_input_cell_array);
%         tmp_input_cell_array = {data_fading_dbm,fading_params,fading_trunc_val_dbm'};
%         ll_fading_dbm_truncated = funonarray(ll_fun_handle,tmp_input_cell_array);
%         % Plot
%         figure;plot(ll_fading_dbm);hold on ;plot(ll_fading_dbm_truncated);grid on;title([minimal_experiment_name,'Log Likelihood RSSI - Distance']);legend('Full Distribution','Truncated Distribution');saveas(gcf,[plot_folder_path,'/','LL RSSI.png']);
%         %% KS-Test
%         h = funoncellarray2input(generated_rssi_dbm_truncated,data_dbm_cell,@(x,y)kstest2(x,y,'Alpha',.01));
%         figure;plot(h);title('KS-Test');xlabel('Distance (m)');ylim([-1,2]);ylabel('Null Hyphotesis State');grid on;saveas(gcf,[plot_folder_path,'/','KS-Test.png']);
end
