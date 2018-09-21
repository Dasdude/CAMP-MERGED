clc
close all
% clear
addpath(genpath('.'));

files_desc = dir('./Dataset/Dataset Raw/**/*.csv');
file_exp_num = zeros(1,length(files_desc));
for i = 1:length(files_desc)
    file_name = files_desc(i).name;
    file_name = split(file_name,'_');
    file_exp_num(i) = str2num(file_name{1});
    file_obe_id(i)= str2num(file_name{2}(2:end));
end
experiments_numbers_set = unique(file_exp_num);
all_merged = table();
for exp_num_idx = 1:length(experiments_numbers_set)
    
    experiment_number = experiments_numbers_set(exp_num_idx);
    experiment_number
    file_masks = file_exp_num==experiment_number;
    same_exp_files = files_desc(file_masks);
    obe_ids = file_obe_id(file_masks);
    all_dataset = cell(1,length(same_exp_files));
    transmitted_dataset= cell(1,length(same_exp_files));
    recieved_dataset= cell(1,length(same_exp_files));
    gps_update_dataset  = cell(1,length(same_exp_files));
    for i = 1:length(same_exp_files)
        sprintf('File %d from %d read',i,length(same_exp_files))
        file_path = fullfile(same_exp_files(i).folder,same_exp_files(i).name);
        all_dataset{i} = {};
        tmp = readtable(file_path,'ReadVariableNames',true);
        all_dataset{i} = tmp(:,{'TimeStamp_ms','lat','long','msgCnt','MsgSeqNum','LogRecType','Range','UniqueOBE_ID_Alias','secMark','RSS'});
        all_dataset{i}.rx_lat = all_dataset{i}.lat;
        all_dataset{i}.rx_long = all_dataset{i}.long;
        all_dataset{i}.rx_lat = fillmissing(all_dataset{i}.lat,'previous');
        all_dataset{i}.rx_long = fillmissing(all_dataset{i}.long,'previous');
        all_dataset{i}.reciever_id = cell2mat(cellfun(@(x)str2num(x(2:end)),all_dataset{i}.UniqueOBE_ID_Alias,'UniformOutput',false));
        tmp = all_dataset{i};
        transmitted_dataset{i} = tmp(strcmp(tmp.LogRecType,'TXE'),:);
        recieved_dataset = tmp(strcmp(tmp.LogRecType,'RXE'),:);
        gps_dataset{i} = tmp(strcmp(tmp.LogRecType,'STE'),:);
    end
    
    all_for_experiment_merged = table();
    for rec_idx = 1:length(same_exp_files)
        
        rec_obe = obe_ids(rec_idx);
        
        for trans_idx = 1:length(same_exp_files)
            
            if trans_idx ==rec_idx
                continue
            end
            sprintf('reciever %d transmitter %d',rec_idx,trans_idx)
            trans_obe = obe_ids(trans_idx);
            rec_table_all = all_dataset{rec_idx};
            trans_table = all_dataset{trans_idx};
            trans_table_test = trans_table(strcmp(trans_table.LogRecType,'TXE'),:);
            trans_table = trans_table_test;
            trans_table{:,'rx_lat'} = nan;
            trans_table{:,'rx_long'} = nan;
            rec_table = rec_table_all(rec_table_all.reciever_id==trans_obe,:);
            rec_table_ste = rec_table_all(strcmp(rec_table_all.LogRecType,'STE'),:);
            rec_table_ste = rec_table_ste(:,{'TimeStamp_ms','lat','long','msgCnt','MsgSeqNum','secMark','Range','rx_lat','rx_long','RSS','LogRecType'});
            trans_table = trans_table(:,{'TimeStamp_ms','lat','long','msgCnt','MsgSeqNum','secMark','Range'});
            tx = trans_table;
            rx = rec_table(:,{'MsgSeqNum','rx_lat','rx_long','RSS','LogRecType','TimeStamp_ms'});             
            joined = outerjoin(tx,rx,'Keys','MsgSeqNum');
            joined_com = sortrows([joined(abs(joined.TimeStamp_ms_rx-joined.TimeStamp_ms_tx)<1000,:);joined(isnan(joined.MsgSeqNum_rx),:)],'TimeStamp_ms_tx');
            joined_com{:,'LogRecType'} = {'RXE'};
            joined_com.MsgSeqNum = joined_com.MsgSeqNum_rx;
            joined_com = removevars(joined_com,{'MsgSeqNum_rx','MsgSeqNum_tx'});
            rec_table_ste{:,'LogRecType'} = {'STE'};
            rec_table_ste.TimeStamp_ms_rx = rec_table_ste.TimeStamp_ms;
            rec_table_ste.TimeStamp_ms_tx = rec_table_ste.TimeStamp_ms;
            rec_table_ste = rec_table_ste(:,joined_com.Properties.VariableNames);
            merged = [rec_table_ste;joined_com];
            merged_sort = sortrows(merged,{'TimeStamp_ms_tx'});
            merged_sort.rx_lat = fillmissing(merged_sort.rx_lat,'previous');
            merged_sort.rx_long = fillmissing(merged_sort.rx_long,'previous');
            merged_recieved = merged_sort(strcmp(merged_sort.LogRecType,'RXE'),:);
            merged_recieved{:,'rx_id'} = rec_obe;
            merged_recieved{:,'tx_id'} = trans_obe;
            d_angle = distance(merged_recieved.rx_lat,merged_recieved.rx_long,merged_recieved.lat,merged_recieved.long);
            range = distdim(d_angle,'deg','meters');
            merged_recieved.Range = range;
            merged_recieved{isnan(merged_recieved.RSS),'RSS'} = -1000;
            merged_recieved.TimeStamp_ms = merged_recieved.TimeStamp_ms_rx;
            merged_recieved = merged_recieved(:,{'TimeStamp_ms','rx_lat','rx_long','lat','long','Range','msgCnt','MsgSeqNum','secMark','RSS','rx_id','tx_id','TimeStamp_ms_rx','TimeStamp_ms_tx'});
            all_for_experiment_merged = [all_for_experiment_merged;merged_recieved];
            all_merged = [all_merged;merged_recieved];
        end
    end
    writetable(all_for_experiment_merged,sprintf('./Dataset/%d_merged.csv',experiment_number));
    
end
writetable(all_merged,sprintf('./Dataset/all_merged.csv'));

%% STEP 2 Separation
data_merged_fixed = clean_intersection_data(all_merged);
writetable(data_merged_fixed,'Dataset/all_merged_rxtxboundry.csv');
[same_leg_flag,nsame_los_flag,nsame_nlos_flag] = insquare(data_merged_fixed);
same_leg_data = data_merged_fixed(same_leg_flag,:);
dif_leg_los_data = data_merged_fixed(nsame_los_flag,:);
dif_leg_nlos_data = data_merged_fixed(nsame_nlos_flag,:);
per_roadmap_heatmap(data_merged_fixed,20,10,'./Plots/DataAnalysis/TRCVR HeatMap/All/');
per_roadmap_heatmap(dif_leg_nlos_data,20,10,'./Plots/DataAnalysis/TRCVR HeatMap/NLOS/');
per_roadmap_heatmap(dif_leg_los,20,10,'./Plots/DataAnalysis/TRCVR HeatMap/LOS/');
per_roadmap_heatmap(same_leg_data,20,10,'./Plots/DataAnalysis/TRCVR HeatMap/SAME/');
writetable(same_leg_data,'Dataset/Seperated/same_leg.csv');
writetable(dif_leg_los_data,'Dataset/Seperated/dif_leg_los.csv');
writetable(dif_leg_nlos_data, 'Dataset/Seperated/dif_leg_nlos.csv');
run('data_analysis_script.m')