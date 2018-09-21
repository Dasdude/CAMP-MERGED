clc
close all

clc
close all
clear
addpath(genpath('.'))
%% File Names
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
for mode_index = 1:3
    mode = mode_list{mode_index};
    mode_name = sprintf('%s',mode{1});
    minimal_experiment_name = [replace(mode{1},'_',' ')];
    plot_folder_path = fullfile('Plots','DataAnalysis','RSS PER',mode{1});
    mkdir(plot_folder_path);
    %% Dataset prepare
    display('Data Prepare Phase')
%     dataset_file_path = sprintf('Dataset/42320000_merged.csv');
    dataset_file_path = sprintf('Dataset/Ehsan/%s.csv',mode{1});
    csv_data = readtable(dataset_file_path,'ReadVariableNames',true);
%     d_max = floor(prctile(csv_data.Range(csv_data.RSS>-300),99.99));
    d_max = 800;
    per = zeros(d_max);
    percentiles_rssi=zeros(d_max,7) ;
    percentiles_rssi_per_inc=zeros(d_max,7) ;
    data_sel = csv_data;
    dataset_mat_dirty = [data_sel.TxRxDistance,data_sel.RSS];
    dataset_mat_dirty(dataset_mat_dirty(:,2)>300,2) = -999;
    dataset_mat_dirty(dataset_mat_dirty(:,2)<-100,2) = -999;
    dataset_cell_dirty = data_mat_cell(dataset_mat_dirty,d_max);
    [dataset_cell,per,packet_loss_stat]=censor_data(dataset_cell_dirty,censor_function_handle);
    packet_received = packet_loss_stat(:,2)-packet_loss_stat(:,1);packet_trans = packet_loss_stat(:,2);
    
    data_dbm_cell = dataset_cell;
    data_dbm_cell = data_dbm_cell(1:d_max);
    data_dbm_mean = funoncellarray1input(data_dbm_cell,@mean);
    data_dbm_std = funoncellarray1input(data_dbm_cell,@std);
    close all     
    %% Percentile Plot
    non_trunc_ylim = [-130,-30];
    percentile_values = [5,10,25,50,75,90,95];
    percentile_values_str =   sprintfc('%d%%',percentile_values);
    percentiles_rssi(:,:) = percentile_array(percentile_values,data_dbm_cell);
    percentiles_rssi_per_inc(:,:) = percentile_array_per([5,10,25,50,75,90,95],data_dbm_cell,per*100);
    figure; plot(squeeze(percentiles_rssi_per_inc(:,:)));grid on;title([minimal_experiment_name,' Percentile']);ylabel('RSS (dbm)');xlabel('Distance (m)');xlim([1,800]);ylim([-100,-30]);legend(percentile_values_str);saveas(gcf,[plot_folder_path,'/','percentile Inc.png']);    
    figure; plot(squeeze(percentiles_rssi(:,:)));grid on;title([minimal_experiment_name,' Percentile']);ylabel('RSS (dbm)');xlabel('Distance (m)');xlim([1,800]);ylim([-100,-30]);legend(percentile_values_str);saveas(gcf,[plot_folder_path,'/','percentile.png']);    
    figure; plot(per*100);grid on;title([minimal_experiment_name,'PER']);ylabel('PER (%)');xlabel('Distance (m)');ylim([0,100]);saveas(gcf,[plot_folder_path,'/','PER.png']);
    figure; plot(medfilt1(per*100,10));grid on;title([minimal_experiment_name,' Smooth PER']);ylabel('PER (%)');ylim([0,100]);xlabel('Distance (m)');saveas(gcf,[plot_folder_path,'/','PER Smooth.png']);
    figure; plot(1:d_max,packet_received,1:d_max,packet_trans);grid on ; title([minimal_experiment_name,' Trans/Rec Samples']);legend('Transmitted','Received');xlabel('Distance (m)');ylabel('# of Samples');saveas(gcf,[plot_folder_path,'/','Samples.png']);
    
end
