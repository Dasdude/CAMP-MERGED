function [data_table_res] = clean_intersection_data(data_table)
%CLEAN_INTERSECTION_DATA Summary of this function goes here
%   Detailed explanation goes here
data_rx_lat = data_table.rx_lat;
data_rx_long = data_table.rx_long;
data_lat = data_table.lat;
data_long = data_table.long;
leftwing_flag = data_rx_lat>33.7585|data_rx_lat<33.7595;
soutwing_flag = data_rx_long>-117.990|data_rx_long<-117.987;
long_range = [-117.990,-117.989];
lat_range = [33.759,33.7595];
% lat_range = [min(data_rx_lat(soutwing_flag)),max(data_rx_lat(soutwing_flag))];
hor_flag = data_lat>lat_range(1) & data_lat<lat_range(2);
ver_flag = data_long>long_range(1)&data_long<long_range(2);
rx_hor_flag = data_rx_lat>lat_range(1) & data_rx_lat<lat_range(2);
rx_ver_flag = data_rx_long>long_range(1)&data_rx_long<long_range(2);

sel_flag = (hor_flag|ver_flag)&(rx_hor_flag|rx_ver_flag);
data_table_res = data_table(sel_flag,:);
figure;
scatter(data_table.lat,data_table.long,1);
title('Transmitter prev')
figure;
scatter(data_table_res.rx_lat,data_table_res.rx_long,1);
title('Receiver')
figure;
scatter(data_table_res.lat,data_table_res.long,1);
title('Transmitter')
end

