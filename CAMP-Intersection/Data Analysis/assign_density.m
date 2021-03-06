clc 
clear
close all
grid on

%% Density Bounds
same_low_up = 15;
same_med_up = 30;
same_high_down = same_med_up;same_low_down = 0;same_high_up = inf;same_med_down = same_low_up;

opposite_low_up = 15;
opposite_med_up = 30;
opposite_high_down = opposite_med_up;opposite_low_down = 0;opposite_high_up = inf;opposite_med_down = opposite_low_up;
%% Data Path
hov_include =0;
addpath(genpath('.'))
addpath(genpath('./..'))
dataset_name = 'DatasetMerged/data_merged.csv';
all_data = readtable(dataset_name,'ReadVariableNames',true);
experiment_name = sprintf('s %d %d o %d %d hov %d',same_low_up,same_med_up,opposite_low_up,opposite_med_up,hov_include);
file_name_string = sprintf('%s/%s/','Seperated DensityPER',experiment_name);
% file_name_string = ['Seperated DensityPER/s';
mkdir(file_name_string)
sanity_check =0;

%% Add Side West is one and East Is zero
all_data.ego_side(:) = 0;
all_data.contact_side(:) = 0;
all_data.ego_side(all_data.rxlat>=33.7743&all_data.rxlong>=-118.06)=1;
all_data.ego_side(all_data.rxlat>=33.77435&all_data.rxlong<=-118.06)=1;
all_data.contact_side(all_data.lat>=33.7743&all_data.long>=-118.06)=1;
all_data.contact_side(all_data.lat>=33.77435&all_data.long<=-118.06)=1;
if hov_include
    all_data.Av_East_Density = all_data.East_Total./8;
    all_data.Av_West_Density = all_data.West_Total./8;
else
    all_data.Av_East_Density = (all_data.x1_east+all_data.x2_east+all_data.x3_east+all_data.x4_east+all_data.x5_east+all_data.x6_east)./6;
    all_data.Av_West_Density = (all_data.x1_west+all_data.x2_west+all_data.x3_west+all_data.x4_west+all_data.x5_west+all_data.x6_west)./6;
end
%% Seperate Encounter Scenario
same_lane_indicator = ((all_data.ego_side==0)&(all_data.contact_side==0))|((all_data.ego_side==1)&(all_data.contact_side==1));
opposite_lane_indicator = ~same_lane_indicator;
same_lane = all_data(same_lane_indicator,:);
opposite_lane = all_data(opposite_lane_indicator,:);

same_lane.average_density(same_lane.ego_side<1)=same_lane.Av_East_Density(same_lane.ego_side<1);
same_lane.average_density(same_lane.ego_side>0)=same_lane.Av_West_Density(same_lane.ego_side>0);
opposite_lane.average_density = (opposite_lane.Av_East_Density+opposite_lane.Av_West_Density)./2;
if sanity_check

    figure;scatter(all_data.rxlat(all_data.rxlat>=33.7743),all_data.rxlong(all_data.rxlat>=33.7743),1,all_data.ego_side(all_data.rxlat>=33.7743));
    figure;scatter(opposite_lane.rxlat,opposite_lane.rxlong,1,opposite_lane.ego_side);
    figure;scatter(opposite_lane.lat,opposite_lane.long,1,opposite_lane.contact_side);
    figure;scatter(same_lane.rxlat,same_lane.rxlong,1,same_lane.ego_side);
    figure;scatter(opposite_lane.lat,opposite_lane.long,1,opposite_lane.contact_side);

end

%% Same Dir Process
same_lane.RSS(same_lane.RSS==-150)=-101;
same_low = same_lane(same_lane.average_density<same_low_up,:);
same_low_clean = same_low(same_low.RSS>-95,:);
same_low_mat = [same_low_clean.Range,same_low_clean.RSS];
same_low_mat_dirty =  [same_low.Range,same_low.RSS];
same_low_cell = data_mat_cell(same_low_mat,800);
same_low_cell_dirty = data_mat_cell(same_low_mat_dirty,800);

same_high = same_lane(same_lane.average_density>=same_high_down,:);
same_high_clean = same_high(same_high.RSS>-95,:);
same_high_mat = [same_high_clean.Range,same_high_clean.RSS];
same_high_mat_dirty = [same_high.Range,same_high.RSS];
same_high_cell = data_mat_cell(same_high_mat,800);
same_high_cell_dirty = data_mat_cell(same_high_mat_dirty,800);

same_med = same_lane(same_lane.average_density>=same_med_down,:);
same_med = same_med(same_med.average_density<same_med_up,:);
same_med_clean = same_med(same_med.RSS>-95,:);
same_med_mat = [same_med_clean.Range,same_med_clean.RSS];
same_med_mat_dirty = [same_med.Range,same_med.RSS];
same_med_cell = data_mat_cell(same_med_mat,800);
same_med_cell_dirty = data_mat_cell(same_med_mat_dirty,800);

same_prc_low = percentile_array([5,20,50,80,95],same_low_cell);
same_prc_low_dirty = percentile_array([5,20,50,80,95],same_low_cell_dirty);
same_prc_low_dirty(same_prc_low_dirty<-100)=nan;

same_prc_med = percentile_array([5,20,50,80,95],same_med_cell);
same_prc_med_dirty = percentile_array([5,20,50,80,95],same_med_cell_dirty);
same_prc_med_dirty(same_prc_med_dirty<-100)=nan;

same_prc_high = percentile_array([5,20,50,80,95],same_high_cell);
same_prc_high_dirty = percentile_array([5,20,50,80,95],same_high_cell_dirty);
same_prc_high_dirty(same_prc_high_dirty<-100)=nan;

same_low_mat_dirty = [same_low.Range,same_low.RSS];
same_low_cell_dirty = data_mat_cell(same_low_mat_dirty,800);
packet_loss_stat_same_low = per_calc(same_low_cell_dirty,-100);
same_per_low = packet_loss_stat_same_low(:,1)./packet_loss_stat_same_low(:,2);
same_total_packets_low = packet_loss_stat_same_low(:,2);
same_recieved_packets_low = packet_loss_stat_same_low(:,2)-packet_loss_stat_same_low(:,1);

same_med_mat_dirty = [same_med.Range,same_med.RSS];
same_med_cell_dirty = data_mat_cell(same_med_mat_dirty,800);
packet_loss_stat_same_med = per_calc(same_med_cell_dirty,-100);
same_total_packets_med = packet_loss_stat_same_med(:,2);
same_recieved_packets_med = packet_loss_stat_same_med(:,2)-packet_loss_stat_same_med(:,1);
same_per_med = packet_loss_stat_same_med(:,1)./packet_loss_stat_same_med(:,2);

same_high_mat_dirty = [same_high.Range,same_high.RSS];
same_high_cell_dirty = data_mat_cell(same_high_mat_dirty,800);
packet_loss_stat_same_high = per_calc(same_high_cell_dirty,-100);
same_total_packets_high = packet_loss_stat_same_high(:,2);
same_recieved_packets_high = packet_loss_stat_same_high(:,2)-packet_loss_stat_same_high(:,1);
same_per_high = packet_loss_stat_same_high(:,1)./packet_loss_stat_same_high(:,2);
%% Opposite Direction
opposite_lane.RSS(opposite_lane.RSS==-150)=-101;
opposite_low = opposite_lane(opposite_lane.average_density<opposite_low_up,:);
opposite_low_clean = opposite_lane(opposite_lane.RSS>-95,:);
opposite_low_mat = [opposite_low_clean.Range,opposite_low_clean.RSS];
opposite_low_mat_dirty =  [opposite_low.Range,opposite_low.RSS];
opposite_low_cell = data_mat_cell(opposite_low_mat,800);
opposite_low_cell_dirty = data_mat_cell(opposite_low_mat_dirty,800);

opposite_med = opposite_lane(opposite_lane.average_density>=opposite_med_down,:);
opposite_med = opposite_med(opposite_med.average_density<opposite_med_up,:);
opposite_med_clean = opposite_med(opposite_med.RSS>-95,:);
opposite_med_mat = [opposite_med_clean.Range,opposite_med_clean.RSS];
opposite_med_mat_dirty = [opposite_med.Range,opposite_med.RSS];
opposite_med_cell = data_mat_cell(opposite_med_mat,800);
opposite_med_cell_dirty = data_mat_cell(opposite_med_mat_dirty,800);

opposite_high = opposite_lane(opposite_lane.average_density>=opposite_high_down,:);
opposite_high_clean = opposite_high(opposite_high.RSS>-95,:);
opposite_high_mat = [opposite_high_clean.Range,opposite_high_clean.RSS];
opposite_high_mat_dirty = [opposite_high.Range,opposite_high.RSS];
opposite_high_cell = data_mat_cell(opposite_high_mat,800);
opposite_high_cell_dirty = data_mat_cell(opposite_high_mat_dirty,800);


opposite_prc_low = percentile_array([5,20,50,80,95],opposite_low_cell);
opposite_prc_low_dirty = percentile_array([5,20,50,80,95],opposite_low_cell_dirty);
opposite_prc_low_dirty(opposite_prc_low_dirty<-100)=nan;

opposite_prc_med = percentile_array([5,20,50,80,95],opposite_med_cell);
opposite_prc_med_dirty = percentile_array([5,20,50,80,95],opposite_med_cell_dirty);
opposite_prc_med_dirty(opposite_prc_med_dirty<-100)=nan;

opposite_prc_high = percentile_array([5,20,50,80,95],opposite_high_cell);
opposite_prc_high_dirty = percentile_array([5,20,50,80,95],opposite_high_cell_dirty);
opposite_prc_high_dirty(opposite_prc_high_dirty<-100)=nan;
%PER
opposite_low_mat_dirty = [opposite_low.Range,opposite_low.RSS];
opposite_low_cell_dirty = data_mat_cell(opposite_low_mat_dirty,800);
packet_loss_stat_opposite_low = per_calc(opposite_low_cell_dirty,-100);
opposite_per_low = packet_loss_stat_opposite_low(:,1)./packet_loss_stat_opposite_low(:,2);
opposite_total_packets_low = packet_loss_stat_opposite_low(:,2);
opposite_recieved_packets_low = packet_loss_stat_opposite_low(:,2)-packet_loss_stat_opposite_low(:,1);

opposite_med_mat_dirty = [opposite_med.Range,opposite_med.RSS];
opposite_med_cell_dirty = data_mat_cell(opposite_med_mat_dirty,800);
packet_loss_stat_opposite_med = per_calc(opposite_med_cell_dirty,-100);
opposite_per_med = packet_loss_stat_opposite_med(:,1)./packet_loss_stat_opposite_med(:,2);
opposite_total_packets_med = packet_loss_stat_opposite_med(:,2);
opposite_recieved_packets_med = packet_loss_stat_opposite_med(:,2)-packet_loss_stat_opposite_med(:,1);

opposite_high_mat_dirty = [opposite_high.Range,opposite_high.RSS];
opposite_high_cell_dirty = data_mat_cell(opposite_high_mat_dirty,800);
packet_loss_stat_opposite_high = per_calc(opposite_high_cell_dirty,-100);
opposite_per_high = packet_loss_stat_opposite_high(:,1)./packet_loss_stat_opposite_high(:,2);
opposite_total_packets_high = packet_loss_stat_opposite_high(:,2);
opposite_recieved_packets_high = packet_loss_stat_opposite_high(:,2)-packet_loss_stat_opposite_high(:,1);
% Density Histogram
figure;histogram(same_high.average_density,0:50);ylim([0,5e5]);xlim([0,50]);hold on;histogram(same_med.average_density,0:50);ylim([0,5e5]);xlim([0,50]);histogram(same_low.average_density,0:50);ylim([0,5e5]);xlim([0,50]);title('Density Histogram for Same Direction');legend('Same High','Same Medium','Same Low');xlabel('Density (Cars/Mile/Lane)');ylabel('Samples');grid on;saveas(gcf,[file_name_string,'/','Same Density Seperate.png']);
figure;histogram(opposite_high.average_density,0:50);ylim([0,5e5]);xlim([0,50]);hold on;histogram(opposite_med.average_density,0:50);ylim([0,5e5]);xlim([0,50]);histogram(opposite_low.average_density,0:50);ylim([0,5e5]);xlim([0,50]);title('Density Histogram for Opposite Direction');legend('Same High','Same Medium','Same Low');xlabel('Density (Cars/Mile/Lane)');ylabel('Samples');grid on;saveas(gcf,[file_name_string,'/','Opposite Density Seperate.png']);
%% Scatter Plot
figure;scatter(same_low_clean.Range,same_low_clean.RSS,.5,same_low_clean.average_density);colorbar();title('Same Low RSSI');xlabel('Range(m)');ylabel('RSSI');grid on;saveas(gcf,[file_name_string,'/','Scatter Same Low RSSI.png']);
figure;scatter(same_med_clean.Range,same_med_clean.RSS,.5,same_med_clean.average_density);colorbar();title('Same Med RSSI');xlabel('Range(m)');ylabel('RSSI');grid on;saveas(gcf,[file_name_string,'/','Scatter Same Med RSSI.png']);
figure;scatter(same_high_clean.Range,same_high_clean.RSS,.5,int64(same_high_clean.average_density)-mod(int64(same_high_clean.average_density),5));colorbar();colormap(parula(5));title('Same High RSSI');xlabel('Range(m)');ylabel('RSSI');grid on;saveas(gcf,[file_name_string,'/','Scatter Same High RSSI.png']);
figure;scatter(opposite_low_clean.Range,opposite_low_clean.RSS,2,opposite_low_clean.average_density);colorbar();title('Opposite Low RSSI');xlabel('Range(m)');ylabel('RSSI');grid on;saveas(gcf,[file_name_string,'/','Scatter Opposite Low RSSI.png']);
figure;scatter(opposite_med_clean.Range,opposite_med_clean.RSS,.5,opposite_med_clean.average_density);colorbar();title('Opposite Med RSSI');xlabel('Range(m)');ylabel('RSSI');grid on;saveas(gcf,[file_name_string,'/','Scatter Opposite Med RSSI.png']);
figure;scatter(opposite_high_clean.Range,opposite_high_clean.RSS,.5,opposite_high_clean.average_density);colorbar();title('Opposite High RSSI');xlabel('Range(m)');ylabel('RSSI');grid on;saveas(gcf,[file_name_string,'/','Scatter Opposite High RSSI.png']);
%% Show Different Percentile
figure;plot(same_prc_low);title('Percentile Same Low Density');xlabel('Range(m)');ylabel('RSSI(dbm)');legend('5%','20%','50%','80%','95%');grid on;saveas(gcf,[file_name_string,'/','Percentile Same Low.png']);
figure;plot(same_prc_low_dirty);title('Percentile Same Low Density Lost Packet Included');xlabel('Range(m)');ylabel('RSSI(dbm)');legend('5%','20%','50%','80%','95%');grid on;saveas(gcf,[file_name_string,'/','Percentile Same Low Lost Packets.png']);
figure;plot(same_prc_med);title('Percentile Same Med Density');xlabel('Range(m)');ylabel('RSSI(dbm)');legend('5%','20%','50%','80%','95%');grid on;saveas(gcf,[file_name_string,'/','Percentile Same Med.png']);
figure;plot(same_prc_med_dirty);title('Percentile Same Med Density Lost Packet Included');xlabel('Range(m)');ylabel('RSSI(dbm)');legend('5%','20%','50%','80%','95%');grid on;saveas(gcf,[file_name_string,'/','Percentile same med Lost Packets.png']);
figure;plot(same_prc_high);title('Percentile Same High Density');legend('5%','20%','50%','80%','95%');xlabel('Range(m)');ylabel('RSSI(dbm)');grid on;saveas(gcf,[file_name_string,'/','Percentile Same High.png']);
figure;plot(same_prc_high_dirty);title('Percentile Same High Density Lost Packet Included');xlabel('Range(m)');ylabel('RSSI(dbm)');legend('5%','20%','50%','80%','95%');grid on;saveas(gcf,[file_name_string,'/','Percentile Same High Lost Packets.png']);

figure;plot(1:800,same_prc_low(:,3),1:800,same_prc_med(:,3),1:800,same_prc_high(:,3));title('Same Low/Med/High Density Comparison');xlabel('Range(m)');ylabel('RSSI(dbm)');legend('Low Density','Medium Density','High Density');grid on;saveas(gcf,[file_name_string,'/','Percentile Same Direction Density Comparison.png']);
figure;plot(1:800,same_prc_low_dirty(:,3),1:800,same_prc_med_dirty(:,3),1:800,same_prc_high_dirty(:,3));title('Same Low/Med/High Density Comparison Lost Packet');xlabel('Range(m)');ylabel('RSSI(dbm)');legend('Low Density','Medium Density','High Density');grid on;saveas(gcf,[file_name_string,'/','Percentile Same Direction Density Comparison Lost Packet.png']);


figure;plot(opposite_prc_low);title('Percentile Opposite Low Density');xlabel('Range(m)');ylabel('RSSI(dbm)');legend('5%','20%','50%','80%','95%');grid on;saveas(gcf,[file_name_string,'/','Percentile opposite Low.png']);
figure;plot(opposite_prc_med);title('Percentile Opposite med Density');xlabel('Range(m)');ylabel('RSSI(dbm)');legend('5%','20%','50%','80%','95%');grid on;saveas(gcf,[file_name_string,'/','Percentile opposite med.png']);
figure;plot(opposite_prc_high);title('Percentile Opposite High Density');xlabel('Range(m)');ylabel('RSSI(dbm)');legend('5%','20%','50%','80%','95%');grid on;saveas(gcf,[file_name_string,'/','Percentile opposite High.png']);
figure;plot(opposite_prc_low_dirty);title('Percentile Opposite Low Density Lost Packet Included');xlabel('Range(m)');ylabel('RSSI(dbm)');legend('5%','20%','50%','80%','95%');grid on;saveas(gcf,[file_name_string,'/','Percentile opposite Low Lost Packets.png']);
figure;plot(opposite_prc_med_dirty);title('Percentile Opposite med Density Lost Packet Included');xlabel('Range(m)');ylabel('RSSI(dbm)');legend('5%','20%','50%','80%','95%');grid on;saveas(gcf,[file_name_string,'/','Percentile opposite med Lost Packets.png']);
figure;plot(opposite_prc_high_dirty);title('Percentile Opposite High Density Lost Packet Included');xlabel('Range(m)');ylabel('RSSI(dbm)');legend('5%','20%','50%','80%','95%');grid on;saveas(gcf,[file_name_string,'/','Percentile opposite High Lost Packets.png']);

figure;plot(1:800,opposite_prc_low(:,3),1:800,opposite_prc_med(:,3),1:800,opposite_prc_high(:,3));xlabel('Range(m)');ylabel('RSSI(dbm)');title('Opposite Low/Med/High Comparison');legend('Low Density','Med Density','High Density');grid on;saveas(gcf,[file_name_string,'/','Percentile Opposite Direction Density Comparison.png']);
figure;plot(1:800,opposite_prc_low_dirty(:,3),1:800,opposite_prc_med_dirty(:,3),1:800,opposite_prc_high_dirty(:,3));xlabel('Range(m)');ylabel('RSSI(dbm)');title('Opposite Low/Med/High Comparison Lost Packet Included');legend('Low Density','Med Density','High Density');grid on;saveas(gcf,[file_name_string,'/','Percentile Opposite Direction Density Comparison Lost Packet.png']);

figure;plot(1:800,same_prc_low_dirty(:,3),1:800,same_prc_high_dirty(:,3),1:800,same_prc_med_dirty(:,3),1:800,opposite_prc_low_dirty(:,3),1:800,opposite_prc_med_dirty(:,3),1:800,opposite_prc_high_dirty(:,3));xlabel('Range(m)');ylabel('RSSI(dbm)');title('Median Same/Op Lo/Med/Hi Density Comparison Lost Packet Included');legend('Same Low','Same Med','Same High','Opposite Low','Opposite Med','Opposite High');grid on;saveas(gcf,[file_name_string,'/','Percentile All Scenario Lost Packet.png']);
%% Show PER Values

figure;plot(1:800,same_per_low,1:800,same_per_high,1:800,same_per_med);title('Same Direction Packet Loss Comparison');xlabel('Range(m)');ylabel('PER');legend('Low Density','High Density','Medium Density','location','northwest');grid on ;saveas(gcf,[file_name_string,'/','PER Same Comparison.png']);
figure;plot(1:800,same_total_packets_high,1:800,same_recieved_packets_high);title('Same High Density Samples Trans/Received');xlabel('Range(m)');ylabel('PER');legend('Total','Received','location','northeast');grid on ;saveas(gcf,[file_name_string,'/','same high samples.png']);
figure;plot(1:800,same_total_packets_low,1:800,same_recieved_packets_low);title('Same Low Density Samples Trans/Received');xlabel('Range(m)');ylabel('PER');legend('Total','Received','location','northeast');grid on ;saveas(gcf,[file_name_string,'/','same low samples.png']);
figure;plot(1:800,same_total_packets_med,1:800,same_recieved_packets_med);title('Same Medium Density Samples Trans/Received');xlabel('Range(m)');ylabel('PER');legend('Total','Received','location','northeast');grid on ;saveas(gcf,[file_name_string,'/','same med samples.png']);

figure;plot(1:800,opposite_per_low,1:800,opposite_per_med,1:800,opposite_per_high);title('Opposite Direction PER Comparison');xlabel('Range(m)');ylabel('PER');legend('Low Density','Medium Density','High Density','location','northwest');grid on ;saveas(gcf,[file_name_string,'/','PER opposite Comparison.png']);
figure;plot(1:800,opposite_total_packets_high,1:800,opposite_recieved_packets_high);title('Opposite High Density Samples Trans/Received');xlabel('Range(m)');ylabel('Samples');legend('Total','Received','location','northeast');grid on ;saveas(gcf,[file_name_string,'/','opposite high samples.png']);
figure;plot(1:800,opposite_total_packets_low,1:800,opposite_recieved_packets_low);title('Opposite Low Density Samples Trans/Received');xlabel('Range(m)');ylabel('Samples');legend('Total','Received','location','northeast');grid on ;saveas(gcf,[file_name_string,'/','opposite low samples.png']);
figure;plot(1:800,opposite_total_packets_med,1:800,opposite_recieved_packets_med);title('Opposite Med Density Samples Trans/Received');xlabel('Range(m)');ylabel('Samples');legend('Total','Received','location','northeast');grid on ;saveas(gcf,[file_name_string,'/','opposite med samples.png']);

figure;plot(1:800,opposite_per_low,1:800,opposite_per_med,1:800,opposite_per_high,1:800,same_per_low,1:800,same_per_med,1:800,same_per_high);title('Same/Opposite PER Comparison');xlabel('Range(m)');ylabel('PER');legend('Op Low Density','Op Medium Density','Op High Density','Same Low Density','Same Medium Density','Same High Density','location','northwest');grid on ;saveas(gcf,[file_name_string,'/','PER all Comparison.png']);
%% Write Data
data_result_folder = sprintf('%s%s',file_name_string,'data results');
mkdir(data_result_folder);
writetable(same_low,[data_result_folder,'/same low.csv']);
writetable(same_high,[data_result_folder,'/same high.csv']);
writetable(same_med,[data_result_folder,'/same med.csv']);

writetable(opposite_low,[data_result_folder,'/opposite low.csv']);
writetable(opposite_high,[data_result_folder,'/opposite high.csv']);
writetable(opposite_med,[data_result_folder,'/opposite med.csv']);