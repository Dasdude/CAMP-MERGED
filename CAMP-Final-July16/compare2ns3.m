clc
close all
clear
addpath(genpath('.'))
experiment_name = 'Camp Highway Model Final July 13 SINR censoring';

noise_level = -98;
pkt_size = -1;
truncation_value = -93;
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
minimal_experiment_name = [upper(mode{1}),' ',upper(mode{2})];
mode_name = [mode{1},' Direction ',mode{2},' Density ',num2str(mode{3}),' to ',num2str(mode{4})];
matlab_plots_folder_path = sprintf('./Plots/%s/%s/',experiment_name,mode_name);
ns3_path = sprintf('./Plots/%s/ns3-channel-model-validation-plots/%s/',experiment_name, mode_name);
ns3_per_osman = get_fig_data(sprintf('%s%s',ns3_path,'aggregatedPER.fig'));
generated_rss_ns3 = get_scatter_fig_data(sprintf('%s%s',ns3_path,'rss_scatter_all.fig'));
censor_function_handle = @(x)censor_function(x,noise_level,pkt_size,truncation_value);
[ns3_gen_trunc,ns3_per] = censor_data(generated_rss_ns3,censor_function_handle);


ns3_gen_trunc_mat = data_cell2mat(ns3_gen_trunc);
ns3_gen_mat = data_cell2mat(generated_rss_ns3);
ns3_per = ns3_per*100;
load(sprintf('%s/%s',matlab_plots_folder_path,'data.mat'));
[generated_rssi_dbm_truncated,generated_per] = censor_data(generated_rssi_dbm,censor_function_handle);
generated_per = generated_per*100;
% tmp_prev_samples = funoncellarray1input(data_dbm_cell,@length);
% total_samples = (100/(100-per)).*tmp_prev_samples;
% 
% 
% [data_dbm_cell,per_add] = censor_data(data_dbm_cell,censor_function_handle);
% tmp_removed_samples = per_add.*tmp_prev_samples;
% per = 100*(1-(tmp_prev_samples-tmp_removed_samples)./total_samples);
%% PER PLOTTING
figure;plot(1:800,generated_per,1:800,ns3_per_osman,1:800,per);grid on;title(sprintf('%s PER Comparison >=%d',minimal_experiment_name,truncation_value));legend('Matlab PER','NS3 PER','Field PER','Location','northwest');saveas(gcf,fullfile(matlab_plots_folder_path,'NS3 PER Sanity Check.png'));
%% PERCENTILE PLOT
percentiles_generated = percentile_array([5,10,25,50,75,90,95],generated_rssi_dbm);
percentiles_generated_trunc = percentile_array([5,10,25,50,75,90,95],generated_rssi_dbm_truncated);
percentiles_rssi = percentile_array([5,10,25,50,75,90,95],data_dbm_cell);
percentiles_gen_ns3 = percentile_array([5,10,25,50,75,90,95],generated_rss_ns3);
percentiles_gen_ns3_trunc = percentile_array([5,10,25,50,75,90,95],ns3_gen_trunc);
% figure;plot(percentiles_generated(:,[1,4,7]));hold on ;plot(percentiles_rssi(:,[1,4,7]));grid on;title([ minimal_experiment_name,' RSS Percentile']);xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');legend('5% model','50% model','95% model','5% field','50% field','95% field');saveas(gcf,['Plots/',file_name_string,'/','Percentile RSS 5.png']);
% figure;plot(percentiles_generated_trunc(:,[1,4,7]));hold on ;plot(percentiles_rssi(:,[1,4,7]));grid on;title([ minimal_experiment_name,' Truncated RSS Percentile']);xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');legend('5% model','50% model','95% model','5% field','50% field','95% field');saveas(gcf,['Plots/',file_name_string,'/','Percentile RSS Truncated 5.png']);
% figure;plot(percentiles_generated(:,[2,4,6]));hold on ;plot(percentiles_rssi(:,[2,4,6]));grid on;title([ minimal_experiment_name,' RSS Percentile']);xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');legend('10% model','50% model','90% model','10% field','50% field','90% field');saveas(gcf,['Plots/',file_name_string,'/','Percentile RSS 10.png']);
% figure;plot(percentiles_generated_trunc(:,[2,4,6]));hold on ;plot(percentiles_rssi(:,[2,4,6]));grid on;title([ minimal_experiment_name,' Truncated RSS Percentile']);xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');legend('10% model','50% model','90% model','10% field','50% field','90% field');saveas(gcf,['Plots/',file_name_string,'/','Percentile RSS Truncated 10.png']);
% figure;plot(percentiles_generated(:,[3,4,5]));hold on ;plot(percentiles_rssi(:,[3,4,5]));grid on;title([ minimal_experiment_name,'  RSS Percentile']);xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');legend('25% model','50% model','75% model','25% field','50% field','75% field');saveas(gcf,['Plots/',file_name_string,'/','Percentile RSS 25.png']);
% figure;plot(percentiles_generated_trunc(:,[3,4,5]));hold on ;plot(percentiles_rssi(:,[3,4,5]));grid on;title([ minimal_experiment_name,' Truncated RSS Percentile']);xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');legend('25% model','50% model','75% model','25% field','50% field','75% field');saveas(gcf,['Plots/',file_name_string,'/','Percentile RSS Truncated 25.png']);
% join('a','b')
figure;plot(percentiles_generated(:,[1,4,7]));hold on ;plot(percentiles_gen_ns3(:,[1,4,7]));grid on;title([ minimal_experiment_name,' RSS Percentile']);xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');legend('5% model','50% model','95% model','5% NS3','50% NS3','95% NS3');saveas(gcf,fullfile(matlab_plots_folder_path,'NS3 Percentile RSS 5.png'));
figure;plot(percentiles_generated_trunc(:,[1,4,7]));hold on ;plot(percentiles_gen_ns3_trunc(:,[1,4,7]));grid on;title([ minimal_experiment_name,' Truncated RSS Percentile']);xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');legend('5% model','50% model','95% model','5% NS3','50% NS3','95% NS3');saveas(gcf,fullfile(matlab_plots_folder_path,'NS3 Percentile RSS Truncated 5.png'));
figure;plot(percentiles_generated(:,[2,4,6]));hold on ;plot(percentiles_gen_ns3(:,[2,4,6]));grid on;title([ minimal_experiment_name,' RSS Percentile']);xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');legend('10% model','50% model','90% model','10% NS3','50% NS3','90% NS3');saveas(gcf,fullfile(matlab_plots_folder_path,'NS3 Percentile RSS 10.png'));
figure;plot(percentiles_generated_trunc(:,[2,4,6]));hold on ;plot(percentiles_gen_ns3_trunc(:,[2,4,6]));grid on;title([ minimal_experiment_name,' Truncated RSS Percentile']);xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');legend('10% model','50% model','90% model','10% NS3','50% NS3','90% NS3');saveas(gcf,fullfile(matlab_plots_folder_path,'NS3 Percentile RSS Truncated 10.png'));
figure;plot(percentiles_generated(:,[3,4,5]));hold on ;plot(percentiles_gen_ns3(:,[3,4,5]));grid on;title([ minimal_experiment_name,'  RSS Percentile']);xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');legend('25% model','50% model','75% model','25% NS3','50% NS3','75% NS3');saveas(gcf,fullfile(matlab_plots_folder_path,'NS3 Percentile RSS 25.png'));
figure;plot(percentiles_generated_trunc(:,[3,4,5]));hold on ;plot(percentiles_gen_ns3_trunc(:,[3,4,5]));grid on;title([ minimal_experiment_name,' Truncated RSS Percentile']);xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');legend('25% model','50% model','75% model','25% NS3','50% NS3','75% NS3');saveas(gcf,fullfile(matlab_plots_folder_path,'NS3 Percentile RSS Truncated 25.png'));

%% Scatter
figure;scatter(ns3_gen_mat(:,1),ns3_gen_mat(:,2),1);grid on;title(sprintf('%s NS3 Scatter Plot',minimal_experiment_name));xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');saveas(gcf,fullfile(matlab_plots_folder_path,'Scatter NS3.png'));
figure;scatter(ns3_gen_trunc_mat(:,1),ns3_gen_trunc_mat(:,2),1);grid on; title(sprintf('%s NS3 Censored Scatter Plot >=%d dbm',minimal_experiment_name,truncation_value));xlabel('TxRx Distance(m)');ylabel('RSS(dbm)');saveas(gcf,fullfile(matlab_plots_folder_path,'Scatter NS3 Censored.png'));
end
function y_data = get_fig_data(fig_path)
% MYMEAN Example of a local function.
    fig_obj = open(sprintf('%s',fig_path));
% fig_obj = open(fig_path);
% close all
    axesObjs = get(fig_obj, 'Children');  %axes handles
    dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
    objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
    xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
    y_data = get(dataObjs, 'YData');
end
function data_cell = get_scatter_fig_data(fig_path);
% MYMEAN Example of a local function.
    fig_obj = open(sprintf('%s',fig_path));
% fig_obj = open(fig_path);
% close all
    axesObjs = get(fig_obj, 'Children');  %axes handles
    dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
%     objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
%     xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
    y_data = dataObjs{3}.YData;
    x_data = dataObjs{3}.XData;
    data_cell = data_mat_cell([x_data',y_data'],800);
end
