function [generated_rssi_dbm] = generated_samples_from_scenario(mode_index,experiment_name,total_samples)
%GENERATED_SAMPLES_FROM_SCENARIO Summary of this function goes here
%   Detailed explanation goes here

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
% experiment_name = 'Camp Highway Model Final July 13 SINR censoring';
minimal_experiment_name = [upper(mode{1}),' ',upper(mode{2})];
mode_name = [mode{1},' Direction ',mode{2},' Density ',num2str(mode{3}),' to ',num2str(mode{4})];
parameter_folder = ['Plots/',experiment_name,'/',mode_name,'/Results'];
parameter_path = [parameter_folder,'/Parameters.mat'];
load(parameter_path);
d_max = 800;
LIGHT_SPEED=3*10^8;
lambda=LIGHT_SPEED/CARRIER_FREQ;
pathloss = pathloss_gen_2ray(TX_HEIGHT,RX_HEIGHT,EPSILON,ALPHA,lambda,d_max);
%% Generate Data

generated_fading_linear = nakagami_generator(fading_params,total_samples);
generated_fading_dbm = linear2dbm(generated_fading_linear);
generated_rssi_dbm = add_fading(pathloss,generated_fading_dbm,TX_POWER);
end

