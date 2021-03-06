clc
close all
clear
addpath(genpath('.'));
dataset_folder = './Dataset';
west_path = [dataset_folder,'/west.csv'];
east_path = [dataset_folder,'/east.csv'];
files = dir('./Dataset/Data-SaifApril/**/*.csv');
total_files = length(files);
%% Crop Range Valuse
min_long = -118.069363;
max_long = -118.051763;
%% Load 
columns_to_store = {'TimeStamp_ms','TimeStamp_S','ego_id','contact_id','rxlat','rxlong','lat','long','RSS','Range','LogRecType','contact_id','TxPwrLevel','speed','heading','LongAccel'};
folder_names = {'March26','March27','March28','March30'};
east = readtable(east_path,'ReadVariableNames',true);
east.Properties.VariableNames{'Total'}='East_Total';
west = readtable(west_path,'ReadVariableNames',true);
west.Properties.VariableNames{'Total'}='West_Total';
density = join(east,west,'keys','TimeStamp_S');
merged = [];
for i = 1:total_files
    file_dir = files(i);
    file_path = [file_dir.folder,'/',file_dir.name];
    
    file_name = file_dir.name;
    veh_id = regexp(file_name,'D0[0-9][0-9]','match');
    veh_id = veh_id{1}(3:end);
    veh_id = str2double(veh_id);
    data_table = readtable(file_path,'ReadVariableNames',true);
    folders_sep = strsplit(file_path,'/');
    
%     data_table = data_table(randi(length(data_table.RSS),20,1),:);
    data_table.ego_id(:) = veh_id;
    data_table = data_table(strcmp(data_table.LogRecType,'RXE'),:);
    data_table.Properties.VariableNames{'UniqueOBE_ID_Alias'} = 'contact_id';
    time_second = int64(data_table.TimeStamp_ms/1000);
    data_table.TimeStamp_S = time_second- mod(time_second,300);
    %% RXLAT AND LONG
    crop_flag = data_table.rxlong>min_long&data_table.rxlong<max_long&data_table.long>min_long&data_table.long<max_long;
    data_crop_ratio = int64(sum(crop_flag)/length(crop_flag)*100);
    fprintf('%s %s : Vehicle %d -  %d in Specified Range - file %d:%d \n',folders_sep{end-2},folders_sep{end-1},veh_id,data_crop_ratio,i,total_files);
    
    data_table_cropped = data_table(crop_flag,:);
    %% Add Density
    data_table_cropped = data_table_cropped(:,columns_to_store);
    q = join(data_table_cropped,density,'keys','TimeStamp_S');
    merged = [merged;q];
%     scatter(merged.rxlat,merged.rxlong)

    
end
merged.ego_side(:) = 0;
merged.contact_side(:) = 0;
merged.ego_side(merged.rxlat>=33.7743&merged.rxlong>=-118.06)=1;
merged.ego_side(merged.rxlat>=33.77435&merged.rxlong<=-118.06)=1;
merged.contact_side(merged.lat>=33.7743&merged.long>=-118.06)=1;
merged.contact_side(merged.lat>=33.77435&merged.long<=-118.06)=1;
% merged.contact_side(merged.lat>=33.774336)=1;
merged.Av_East_Density = merged.East_Total./8;
merged.Av_West_Density = merged.West_Total./8;
%% Seperate Encounter Scenario
same_lane_indicator = ((merged.ego_side==0)&(merged.contact_side==0))|((merged.ego_side==1)&(merged.contact_side==1));
opposite_lane_indicator = ~same_lane_indicator;
same_lane = merged(same_lane_indicator,:);
opposite_lane = merged(opposite_lane_indicator,:);

same_lane.average_density(same_lane.ego_side<1)=same_lane.Av_East_Density(same_lane.ego_side<1);
same_lane.average_density(same_lane.ego_side>0)=same_lane.Av_West_Density(same_lane.ego_side>0);
opposite_lane.average_density = (opposite_lane.Av_East_Density+opposite_lane.Av_West_Density)./2;
same_lane.same_lane_indicator(:) = 1;
opposite_lane.same_lane_indicator(:)=0;
file_to_store_name = ['processed/data_merged.csv'];
writetable(merged,file_to_store_name);