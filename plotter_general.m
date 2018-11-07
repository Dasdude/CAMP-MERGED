clc
close all

read_file_flag = 1;
addpath(genpath('.'))

set(groot,'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0, 'DefaultFigureVisible', 'off');
set(0, 'DefaultLineLineWidth', .5);
%% File Names
PLOT_EDITED_PARAMS = 0;
% dataset_name = 'Intersection';
dataset_folders = dir('Dataset');
for file_idx =3:length(dataset_folders)
    display(sprintf('[%d] %s',file_idx-2,dataset_folders(file_idx).name))
end
dataset_idx = str2num(input('dataset idx?','s'))+2;
clc
dataset_name = dataset_folders(dataset_idx).name;
plot_folders = dir(fullfile('Plots',dataset_name));
for file_idx =3:length(plot_folders)
    display(sprintf('[%d] %s',file_idx-2,plot_folders(file_idx).name))
end
experiment_idx = str2num(input('Experiment idx?','s'))+2;
if isempty(experiment_idx)
    experiment_name = 'Debug';
else
    experiment_name = plot_folders(experiment_idx).name;
end
clc

dataset_folder_path = fullfile('./Dataset',dataset_name);
dataset_files = dir(fullfile(dataset_folder_path,'**/*.csv'));
for file_idx =1:length(dataset_files)
    file = dataset_files(file_idx);
    [~,dataset_files(file_idx).name_wo_extension,~] = fileparts(file.name);
    dataset_files(file_idx).path = fullfile(file.folder,file.name);
end

addpath(genpath('.'))
%% Distribution
dist_obj_nak = distribution_type_class(@(x)makedist('nakagami',x(1),x(2)),'nakagami',{'mu','omega'},[0.5,0],[inf,inf]);
dist_obj_logn = distribution_type_class(@(x)makedist('lognormal',x(1),x(2)),'lognormal',{'mu','sigma'},[-inf,0],[3,3]);
dist_obj_wei = distribution_type_class(@(x)makedist('weibull',x(1),x(2)),'weibull',{'A','B'},[0,0],[5,5]);
dist_obj_ri = distribution_type_class(@(x)makedist('rician',x(1),x(2)),'rician',{'s','sigma'},[0,0],[5,5]);
dist_obj_ray = distribution_type_class(@(x)makedist('rayleigh',x(1)),'rayleigh',{'B'},[0],[5]);
dist_obj_cell = {dist_obj_logn,dist_obj_wei,dist_obj_nak};
dist_cellname = {'Log-Normal','Weibull','Nakagami'};

if PLOT_EDITED_PARAMS
    output_extension = '_edit.png';
    param_file_name = 'Parameters_edit.mat';
else
    output_extension_list = {'','.fig','.png'};
    param_file_name = 'Parameters.mat';
end
noise_level = -98;
pkt_size = -1;
TRUNCATION_VALUE = -90;
med_filt_size = 10;
censor_function_handle = @(x)censor_function(x,noise_level,pkt_size,TRUNCATION_VALUE);
for mode_index = 1:length(dataset_files)
%     mode = mode_list{mode_index};
    read_file_flag =1;
    
    data_file_obj = dataset_files(mode_index);
    %% Dataset prepare
    sprintf('Data Prepare Phase')
    
    if logical(exist('dataset_file_path'))&&strcmp(dataset_file_path,data_file_obj.path)
        read_file_flag = 0;
    end
    dataset_file_path  = data_file_obj.path;
    csv_data = readtable(dataset_file_path,'ReadVariableNames',true);
    if sum(strcmp(csv_data.Properties.VariableNames,'TxRxDistance'))
            csv_data.Range = csv_data.TxRxDistance;
            csv_data.TxRxDistance = [];
            writetable(csv_data,dataset_file_path);
    end
    csv_data.RSS(csv_data.RSS==999) = -999;
    ll = zeros(2000,length(dist_obj_cell));
    ll_dscrt = zeros(2000,length(dist_obj_cell));
    ll_rate = zeros(2000,length(dist_obj_cell));
    model_per = zeros(2000,length(dist_obj_cell));
    for dist_index = 1:length(dist_cellname)      
        close all
        dist_obj = dist_obj_cell{dist_index};
        dist_name = dist_obj.dist_name;dist_name = [upper(dist_name(1)),dist_name(2:end)];
        minimal_experiment_name = [data_file_obj.name_wo_extension,' ',dist_name,' '];
        sprintf('%s',minimal_experiment_name)
        minimal_experiment_name_all_dist = [data_file_obj.name_wo_extension,' '];
        relative_experiment_folder_path = fullfile(dataset_name,experiment_name,data_file_obj.name_wo_extension,dist_obj.dist_name);
        parameter_folder = fullfile(relative_experiment_folder_path,'Results');
        parameter_path = fullfile(parameter_folder,param_file_name);
        plot_folder_path = fullfile('./Plots',relative_experiment_folder_path);
        %% Load Params
        
        try
            load(parameter_path);
        catch err
            display(err)
            continue
        end
        fading_params = medfilt1(fading_params,med_filt_size);
        ll = ll(1:size(fading_params,1),:);
        ll_dscrt = ll_dscrt(1:size(fading_params,1),:);
        ll_rate = ll_rate(1:size(fading_params,1),:);
        model_per = model_per(1:size(fading_params,1),:);
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
        data_fading_linear = dbm2linear(data_fading_dbm);
        
        fading_truncation_val_dbm =  TRUNCATION_VALUE-TX_POWER+pathloss;
        fading_trunc_val_linear = dbm2linear(fading_truncation_val_dbm);
        %% Generate Data
        generated_fading_linear = sample_generator(dist_obj,fading_params,1e3);
        generated_fading_dbm = linear2dbm(generated_fading_linear);
        generated_rssi_dbm = add_fading(pathloss,generated_fading_dbm,TX_POWER);
        [generated_rssi_dbm_truncated,generated_per,gen_pl_stat] = censor_data(generated_rssi_dbm,censor_function_handle);
        generated_rssi_dbm_mean = funoncellarray1input(generated_rssi_dbm,@mean);
        model_per(:,dist_index) = generated_per';
           %% Pathloss Compare Plot
        figure;subplot(2,1,1);plot(generated_rssi_dbm_mean);hold;plot(data_dbm_mean);title([minimal_experiment_name,'Mean Comparison']);grid on;legend('Model','Field');subplot(2,1,2);plot(aprx_per);title('PER');ylabel('RSS');saveas(gcf,[plot_folder_path,'/','Mean Model vs Field',output_extension_list{1}]);
        %% Percentile Plot
        non_trunc_ylim = [-130,-30];
        percentile_generation_values = [5,10,25,50,75,90,95];
        percentile_string_array = sprintfc('%d %%',percentile_generation_values);
        percentiles_generated = percentile_array(percentile_generation_values,generated_rssi_dbm);
        percentiles_generated_trunc = percentile_array(percentile_generation_values,generated_rssi_dbm_truncated);
        percentiles_rssi = percentile_array(percentile_generation_values,data_dbm_cell);
        percentiles_rssi_per_inc = percentile_array_per(percentile_generation_values,data_dbm_cell,per*100);
        percentiles_rssi_gen_per_inc = percentile_array_per(percentile_generation_values,generated_rssi_dbm_truncated,generated_per*100);
        save([plot_folder_path,'/','Datamat.mat'],'percentiles_generated','percentiles_generated_trunc','percentiles_rssi','percentiles_rssi_per_inc','percentiles_rssi_gen_per_inc','generated_per','per','percentile_generation_values','pathloss');
        for oe_idx = 1:length(output_extension_list)
        output_extension = output_extension_list{oe_idx};
%         indices = [2,4,6];percentile_plot(indices,'Percentile 10 Full',percentiles_generated,percentiles_rssi_per_inc,percentile_generation_values,minimal_experiment_name,non_trunc_ylim,plot_folder_path,output_extension);
%         indices = [3,4,5];percentile_plot(indices,'Percentile 25 Full',percentiles_generated,percentiles_rssi_per_inc,percentile_generation_values,minimal_experiment_name,non_trunc_ylim,plot_folder_path,output_extension);
%         indices = [1,4,7];percentile_plot(indices,'Percentile 5 Full',percentiles_generated,percentiles_rssi_per_inc,percentile_generation_values,minimal_experiment_name,non_trunc_ylim,plot_folder_path,output_extension);
        indices = [1,3,4,5,7];percentile_plot(indices,'Percentile ALL Full',percentiles_generated,percentiles_rssi_per_inc,percentile_generation_values,minimal_experiment_name,non_trunc_ylim,plot_folder_path,output_extension);
        indices = [1,3,4,5,7];percentile_plot(indices,'Percentile ALL Full Trunc',percentiles_generated_trunc,percentiles_rssi,percentile_generation_values,minimal_experiment_name,non_trunc_ylim,plot_folder_path,output_extension);
        %% PER Plot
        figure;plot(packet_loss_stat(:,2));hold;plot(packet_loss_stat(:,2)-packet_loss_stat(:,1));xlabel('Distance(m)');ylabel('Number of Samples');grid on;title([minimal_experiment_name,'Total Samples vs Received Samples']);legend('Total Samples','Recieved Samples');saveas(gcf,[plot_folder_path,'/','Samples Received vs Total',output_extension]);
        figure; plot(100*generated_per);hold on;plot(100*per);grid on;title([minimal_experiment_name,'PER Value Comparison']);ylabel('Rate');xlabel('Distance (m)');legend('Model','Field','Location','best');saveas(gcf,[plot_folder_path,'/','PER Comparison',output_extension]);
        figure;plot(loss_vals);title([minimal_experiment_name,'loss']);grid on;saveas(gcf,[plot_folder_path,'/','Loss',output_extension]);
        %% Plot Parameters
        for param_idx = 1:dist_obj.get_dof
            param_name = dist_obj.dist_params_names{param_idx};
            figure;plot(fading_params(:,param_idx));title(sprintf('%s Parameter: %s - Distance',minimal_experiment_name,param_name));grid on;xlabel('Distance (m)');ylabel(sprintf('%s Value',param_name));saveas(gcf,fullfile(plot_folder_path,sprintf('%s_distance%s',param_name,output_extension)));
        end
        %% LOG Likelihood
        [ll_temp,ll_d,ll_rate_temp]=general_loglikelihood(dist_obj,fading_params,data_fading_linear,fading_trunc_val_linear);
%         ll_temp = medfilt1(ll_temp,5,'omitnan');
%         ll_d = medfilt1(ll_d,5,'omitnan');
        ll_rate_temp = medfilt1(ll_rate_temp,5,'omitnan');
        figure;plot(ll_temp);title([minimal_experiment_name,' Log Likelihood Descreet']);xlabel('Distance (m)');ylabel('Average LL');grid on;saveas(gcf,[plot_folder_path,'/','Log Likelihood Descreet',output_extension]);
        figure;plot(ll_d);title([minimal_experiment_name,' Log Likelihood']);grid on;xlabel('Distance (m)');ylabel('Average LL');saveas(gcf,[plot_folder_path,'/','Log Likelihood',output_extension]);
        figure;plot(ll_d);title([minimal_experiment_name,' Log Likelihood Rate']);grid on;saveas(gcf,[plot_folder_path,'/','Log Likelihood Rate',output_extension]);
        ll(:,dist_index)=ll_temp';
        ll_dscrt(:,dist_index)=ll_d';
        ll_rate(:,dist_index) = ll_rate_temp';
        
%         funoncellarray1input(generated_rssi_dbm_truncated,@floor);
%         funoncellarray2input(generated_rssi_dbm_truncated,   data_dbm_cell,@kstest2)
        end
    end
    for oe_idx = 1:length(output_extension_list)
        output_extension = output_extension_list{oe_idx};
        plot_folder_path_all_dist = fullfile('./Plots',dataset_name,experiment_name,data_file_obj.name_wo_extension);
        ll_dscrt(ll_dscrt>10)=nan;
        figure;plot(ll);title([minimal_experiment_name_all_dist,' Log Likelihood']);xlabel('Distance(m)');ylabel('Average Log Likelihood');legend(dist_cellname);grid on;saveas(gcf,[plot_folder_path_all_dist,'/','ALL Log Likelihood',minimal_experiment_name_all_dist,output_extension]);
        figure;plot(packet_loss_stat(:,2));hold;plot(packet_loss_stat(:,2)-packet_loss_stat(:,1));xlabel('Distance(m)');ylabel('Number of Samples');grid on;title([minimal_experiment_name_all_dist,'Total Samples vs Received Samples']);legend('Total Samples','Recieved Samples');saveas(gcf,[plot_folder_path_all_dist,'/','Samples Received vs Total ALL',output_extension]);
        figure;plot(ll_dscrt);title([minimal_experiment_name_all_dist,' Log Likelihood Discrete']);ylim([0,10]);legend(dist_cellname,'Location','best');grid on;xlabel('Distance(m)');ylabel('Average Log Likelihood');saveas(gcf,[plot_folder_path_all_dist,'/','ALL Log Likelihood Descreet',minimal_experiment_name_all_dist,output_extension]);
        figure;plot(medfilt1(ll_dscrt,5,'omitnan'));title([minimal_experiment_name_all_dist,' Log Likelihood Discrete']);ylim([0,10]);legend(dist_cellname,'Location','best');grid on;xlabel('Distance(m)');ylabel('Average Log Likelihood');saveas(gcf,[plot_folder_path_all_dist,'/','ALL Smooth Log Likelihood Descreet',minimal_experiment_name_all_dist,output_extension]);
        per_labels = dist_cellname;per_labels(end+1) = {'Field'};
        figure;plot([100*model_per,100*per']);title([minimal_experiment_name_all_dist,' PER']);legend(per_labels,'Location','best');grid on;xlabel('Distance(m)');ylabel('PER Rate (\%)');saveas(gcf,[plot_folder_path_all_dist,'/','ALL PER Models',minimal_experiment_name_all_dist,output_extension]);
        figure;plot(medfilt1([100*model_per,100*per'],5,'omitnan'));title([minimal_experiment_name_all_dist,' PER']);legend(per_labels,'Location','best');grid on;xlabel('Distance(m)');ylabel('PER Rate (\%)');saveas(gcf,[plot_folder_path_all_dist,'/','ALL Smooth PER Models',minimal_experiment_name_all_dist,output_extension]);
    end
end
function [] = percentile_plot(indices,file_name,percentile_model,percentile_field,percentile_generation_values,minimal_experiment_name,non_trunc_ylim,plot_folder_path,output_extension)
        default_blue = [0, 0.4470, 0.7410];
        default_orange = [0.8500, 0.3250, 0.0980];
        percentile_string_array = sprintfc('%d%%',percentile_generation_values(indices));
        str_perc = strjoin(percentile_string_array);
        str_model =  strcat(str_perc,'Model');
        str_field =  strcat(str_perc,'Field');
        pstr_model = strcat(percentile_string_array,' Model');
        pstr_field = strcat(percentile_string_array,' Field');
        percentiles_generated = percentile_model;
        percentiles_rssi_per_inc = percentile_field;
        figure;plot_handle_m = plot(percentiles_generated(:,indices),'Color',default_blue);hold on ;plot_handle_f = plot(percentiles_rssi_per_inc(:,indices),'Color',default_orange);
        grid on;
        ylim(non_trunc_ylim);
        title([minimal_experiment_name,'Percentile']);
        ylabel('RSS (dbm)','interpreter','latex');
        xlabel('Distance (m)');
        h = zeros(2, 1);
        h(1) = plot(NaN,NaN,'Color',default_blue);
        h(2) = plot(NaN,NaN,'Color',default_orange);

        legend(h, str_model,str_field,'Interpreter','latex');
%         legend(pstr_model,pstr_field)
%         legend([plot_handle_m,plot_handle_f],{'Model','Field'});
        saveas(gcf,[plot_folder_path,'/',file_name,minimal_experiment_name,  output_extension]);
        
end