clc
close all
% clear
addpath(genpath('.'));
target_density = 15;
total_samples = 1e3;
censor_noise_level = -98;
censor_pkt_size = -1;
censor_trunc_val = -93;
opposite_dir_flag =0;
a_sigmoid = .3;
bias_sigmoid = 22;
%% 
experiment_name = 'Camp Highway Model Final July 30- 14dbm';
if opposite_dir_flag
    scenario ='opposite';
else
    scenario='same';
end
plot_dir_folder =fullfile('Plots',experiment_name,'interpolation Params',scenario);
mkdir(plot_dir_folder)
%%
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
mode = mode_list{3*opposite_dir_flag+1};
mode_name = [mode{1},' Direction ',mode{2},' Density ',num2str(mode{3}),' to ',num2str(mode{4})];
parameter_folder = ['Plots/',experiment_name,'/',mode_name,'/Results'];
parameter_path = [parameter_folder,'/Parameters.mat'];
load(parameter_path);
low_fading_params = fading_params;
mode = mode_list{3*opposite_dir_flag+3};
mode_name = [mode{1},' Direction ',mode{2},' Density ',num2str(mode{3}),' to ',num2str(mode{4})];
parameter_folder = ['Plots/',experiment_name,'/',mode_name,'/Results'];
parameter_path = [parameter_folder,'/Parameters.mat'];
load(parameter_path);
high_fading_params = fading_params;
%%


same_low_rss = generated_samples_from_scenario(1+(3*opposite_dir_flag),'Camp Highway Model Final July 13 SINR censoring',int64(total_samples));
same_high_rss = generated_samples_from_scenario(3+(3*opposite_dir_flag),'Camp Highway Model Final July 13 SINR censoring',int64(total_samples));
same_high_rss_mat = data_cell2mat(same_high_rss);
same_low_rss_mat = data_cell2mat(same_low_rss);
field_med_dirty = load_field_data_scenario(scenario,'med',800);
censor_function_handle = @(x)censor_function(x,censor_noise_level,censor_pkt_size,censor_trunc_val);
[field_med_dataset_cell,field_med_per] = censor_data(field_med_dirty,censor_function_handle);
field_med_prctl = percentile_array([10,50,90],field_med_dataset_cell);
param_fig = figure;
rss_fig = figure('Position',[1 1 2000 1000]);
total_plots = 6;
densities = linspace(0,50,total_plots);
densities = [10,25,40];
total_plots = length(densities);
for i = 1:total_plots
    target_density = densities(i);
    p_high = sigmoid(target_density,a_sigmoid,bias_sigmoid);
    interp_fading_params = p_high*high_fading_params+((1-p_high)*low_fading_params);
    generated_med_cell =generated_samples_from_parameters(TX_HEIGHT,RX_HEIGHT,EPSILON,CARRIER_FREQ,ALPHA,interp_fading_params,TX_POWER,800,total_samples);
%     generated_med_cell = mix_data_cell_array(same_high_rss,same_low_rss,p_high);
    [generated_med_cell_trunc,gen_per] = censor_data(generated_med_cell,censor_function_handle);
    gen_trunc = percentile_array([10,50,90],generated_med_cell_trunc);
    figure(rss_fig);
    subplot(2,total_plots,i);plot(1:800,gen_trunc,1:800,field_med_prctl);grid on;legend('10% gen','50% gen','90% gen','10% field','50% field','90% field');title(sprintf('d:%g p:%g',target_density,p_high));
    subplot(2,total_plots,i+total_plots);plot(1:800,field_med_per,1:800,gen_per);grid on;legend('Field Med Density PER','Generated PER');ylabel('PER');xlabel('TxRxRange(m)');
    figure(param_fig);
    
    subplot(2,total_plots,i);plot(1:800,interp_fading_params(:,1));grid on;
    subplot(2,total_plots,i+total_plots);plot(1:800,interp_fading_params(:,2));grid on;
    
    
end
suptitle(sprintf('%s Param Interpolation a:%g b:%g',upper(scenario), a_sigmoid,bias_sigmoid));
saveas(rss_fig,fullfile(plot_dir_folder,sprintf('rss %d.jpg',total_plots)));
saveas(param_fig,fullfile(plot_dir_folder,sprintf('param %d.jpg',total_plots)));

