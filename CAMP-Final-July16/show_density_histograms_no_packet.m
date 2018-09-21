clc
close all
clear
include_packet_loss = 1;
if include_packet_loss
    folder_path = fullfile('Plots','Data Analysis','Veh Density Lost Packets');
else
    folder_path = fullfile('Plots','Data Analysis','Veh Density Without Lost Packets');
end

mkdir(folder_path);
if ~exist('dens_east')
    [field_same_low_dirty,field_same_low_dirty_mat,dens_east_same_low,dens_west_same_low,dens_same_low] = load_field_data_scenario('same','low',800,include_packet_loss);
    [field_same_med_dirty,field_same_med_dirty_mat,dens_east_same_med,dens_west_same_med,dens_same_med] = load_field_data_scenario('same','med',800,include_packet_loss);
    [field_same_high_dirty,field_same_high_dirty_mat,dens_east_same_high,dens_west_same_high,dens_same_high] = load_field_data_scenario('same','high',800,include_packet_loss);

    [field_opposite_low_dirty,field_opposite_low_dirty_mat,dens_east_op_low,dens_west_op_low,~] = load_field_data_scenario('opposite','low',800,include_packet_loss);
    [field_opposite_med_dirty,field_opposite_med_dirty_mat,dens_east_op_med,dens_west_op_med,~] = load_field_data_scenario('opposite','med',800,include_packet_loss);
    [field_opposite_high_dirty,field_opposite_high_dirty_mat,dens_east_op_high,dens_west_op_high,~] = load_field_data_scenario('opposite','high',800,include_packet_loss);
    dens_east = [dens_east_same_high;dens_east_same_med;dens_east_same_low;dens_east_op_low;dens_east_op_med;dens_east_op_high];
    dens_west = [dens_west_same_high;dens_west_same_med;dens_west_same_low;dens_west_op_low;dens_west_op_med;dens_west_op_high];
    dens_op_west = [dens_west_op_low;dens_west_op_med;dens_west_op_high];
    dens_op_east = [dens_east_op_low;dens_east_op_med;dens_east_op_high];
    dens_same_west = [dens_west_same_low;dens_west_same_med;dens_west_same_high];
    dens_same_east = [dens_east_same_low;dens_east_same_med;dens_east_same_high];
    dens_same = [dens_same_low;dens_same_med;dens_same_high];
    dens_op = [dens_op_west,dens_op_east];
    dens_op_low = [dens_east_op_low,dens_west_op_low];
    dens_op_med = [dens_east_op_med,dens_west_op_med];
    dens_op_high = [dens_east_op_high,dens_west_op_high];
%     dens_same = [dens_same_west,dens_same_east];
    dens_same_no_hov = dens_same(:,1:6);
    dens_op_no_hov = [dens_op_west(:,1:6),dens_op_east(:,1:6)];
end
% plot_lane_wise(dens_same,folder_path,'Following');
% plot_lane_wise(dens_east,folder_path,'North');
% plot_lane_wise(dens_west,folder_path,'South');
figure('Position',[1,1,2000,1000]);histogram(table2array([dens_east,dens_west]));title('All Agg Veh Density (V/M/L)');saveas(gcf,fullfile(folder_path,'AllAgg.png'));
figure('Position',[1,1,2000,1000]);histogram(table2array(dens_op));title('Encounter Scenario Agg Veh Density (V/M/L)');saveas(gcf,fullfile(folder_path,'HOV EncounerAgg.png')); 
figure('Position',[1,1,2000,1000]);histogram(table2array(dens_same));title('Following Scenario Agg Veh Density (V/M/L)');saveas(gcf,fullfile(folder_path,'HOV FollowingAgg.png')); 
figure('Position',[1,1,2000,1000]);histogram(table2array(dens_same_no_hov));title('Following Scenario Agg Veh Density (V/M/L) No HOV');saveas(gcf,fullfile(folder_path,'No HOV FollwoingAgg.png')); 
figure('Position',[1,1,2000,1000]);histogram(table2array(dens_op_no_hov));title('Encounter Scenario Agg Veh Density (V/M/L) No HOV');saveas(gcf,fullfile(folder_path,'NO HOV EncounerAgg.png')); 

figure('Position',[1,1,2000,1000]);histogram(mean(table2array(dens_op),2));title('Encounter Scenario Average Veh Density (V/M/L)');saveas(gcf,fullfile(folder_path,'HOV EncounerAverage.png')); 
figure('Position',[1,1,2000,1000]);histogram(mean(table2array(dens_same),2));title('Following Scenario Average Veh Density (V/M/L)');saveas(gcf,fullfile(folder_path,'HOV FollowingAverage.png')); 
figure('Position',[1,1,2000,1000]);histogram(mean(table2array(dens_same_no_hov),2));title('Following Scenario Average Veh Density (V/M/L) No HOV');saveas(gcf,fullfile(folder_path,'No HOV FollwoingAverage.png')); 
figure('Position',[1,1,2000,1000]);histogram(mean(table2array(dens_op_no_hov),2));title('Encounter Scenario Average Veh Density (V/M/L) No HOV');saveas(gcf,fullfile(folder_path,'NO HOV EncounerAverage.png')); 


% figure('Position',[1,1,2000,1000]);histogram(mean(table2array([dense_east_op_low(:,1:6),dense_west_op_low(:,1:6)]),2));title('Encounter Low Density Scenario Average Veh Density (V/M/L) No HOV');saveas(gcf,fullfile(folder_path,'NO HOV EncounerAverage Low.png')); 
% figure('Position',[1,1,2000,1000]);histogram(mean(table2array([dense_east_op_med(:,1:6),dense_west_op_med(:,1:6)]),2));title('Encounter Med Density Scenario Average Veh Density (V/M/L) No HOV');saveas(gcf,fullfile(folder_path,'NO HOV EncounerAverage med.png')); 
% figure('Position',[1,1,2000,1000]);histogram(mean(table2array([dense_east_op_high(:,1:6),dense_west_op_high(:,1:6)]),2));title('Encounter High Density Scenario Average Veh Density (V/M/L) No HOV');saveas(gcf,fullfile(folder_path,'NO HOV EncounerAverage high.png')); 
%% Cluster Density Same No _HOV
figure('Position',[1,1,2000,1000]);histogram(mean(table2array(dens_same_low(:,1:6)),2),'BinEdges',[0:60]);title('Following Scenario Average Veh Density (V/M/L) No HOV');ylim([0 3e5]);xlim([0,60]);
hold on;histogram(mean(table2array(dens_same_med(:,1:6)),2),'BinEdges',[0:60]);ylim([0 3e5]);xlim([0,60]);
hold on;histogram(mean(table2array(dens_same_high(:,1:6)),2),'BinEdges',[0:60]);ylim([0 3e5]);xlim([0,60]);legend('Low','Medium','High');saveas(gcf,fullfile(folder_path,'NO HOV FollowingAverage Cluster.png')); 
%% Cluster Density Opposite NO HOV

figure('Position',[1,1,2000,1000]);histogram(mean(table2array(dens_op_low(:,[1:6,9:end])),2),'BinEdges',[0:60]);title('Encounter  Scenario Average Veh Density (V/M/L) No HOV');ylim([0 3e5]);xlim([0,50]);
hold on;histogram(mean(table2array(dens_op_med(:,[1:6,9:end])),2),'BinEdges',[0:60]);ylim([0 3e5]);xlim([0,60]);
hold on;histogram(mean(table2array(dens_op_high(:,[1:6,9:end])),2),'BinEdges',[0:60]);ylim([0 3e5]);xlim([0,60]);legend('Low','Medium','High');saveas(gcf,fullfile(folder_path,'NO HOV EncounterAverage Cluster.png')); 
%% Cluster Density Same  HOV

figure('Position',[1,1,2000,1000]);histogram(mean(table2array(dens_same_low),2),'BinEdges',[0:60]);title('Following Scenario Average Veh Density (V/M/L) HOV Included');ylim([0 3e5]);xlim([0,50]);
hold on;histogram(mean(table2array(dens_same_med),2),'BinEdges',[0:60]);ylim([0 3e5]);xlim([0,60]);
hold on;histogram(mean(table2array(dens_same_high),2),'BinEdges',[0:60]);ylim([0 3e5]);xlim([0,60]);legend('Low','Medium','High');saveas(gcf,fullfile(folder_path,'HOV Included FollowingAverage Cluster.png')); 
%% Cluster Density Opposite HOV
figure('Position',[1,1,2000,1000]);histogram(mean(table2array(dens_op_low),2),'BinEdges',[0:60]);title('Encounter Scenario Average Veh Density(V/M/L) HOV Included');ylim([0 3e5]);xlim([0,50]);
hold on;histogram(mean(table2array(dens_op_med),2),'BinEdges',[0:60]);ylim([0 3e5]);xlim([0,60]);
hold on;histogram(mean(table2array(dens_op_high),2),'BinEdges',[0:60]);ylim([0 3e5]);xlim([0,60]);legend('Low','Medium','High');saveas(gcf,fullfile(folder_path,'HOV Included EncounterAverage Cluster.png')); 


function [] = plot_lane_wise(dens,folder_path,bound_name)
    figure('Position',[1,1,4000,4000]);
    values = dens.Variables;
    for i = 1:length(dens.Properties.VariableNames)
        subplot(3,3,i)

        histogram(values(:,i));
        title(sprintf('Lane %g',i));

    end
    subplot(3,3,9);
    histogram(values(:));
    title(sprintf('All'));
    suptitle(sprintf('%s Bound Densities(Veh/Mile/Lane)',bound_name));
    saveas(gcf,fullfile(folder_path,[bound_name,'.png']))
end
