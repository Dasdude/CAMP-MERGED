function [generated_rssi_dbm] = generated_samples_from_parameters(TX_HEIGHT,RX_HEIGHT,EPSILON,CARRIER_FREQ,ALPHA,fading_params,TX_POWER,d_max,total_samples)
%GENERATED_SAMPLES_FROM_PARAMETERS Summary of this function goes here
%   Detailed explanation goes here
LIGHT_SPEED=3*10^8;
lambda=LIGHT_SPEED/CARRIER_FREQ;
pathloss = pathloss_gen_2ray(TX_HEIGHT,RX_HEIGHT,EPSILON,ALPHA,lambda,d_max);
%% Generate Data

generated_fading_linear = nakagami_generator(fading_params,total_samples);
generated_fading_dbm = linear2dbm(generated_fading_linear);
generated_rssi_dbm = add_fading(pathloss,generated_fading_dbm,TX_POWER);
end

