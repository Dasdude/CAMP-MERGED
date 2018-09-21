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
plot_dir_folder =fullfile('Plots',experiment_name,'interpolation Mixture',scenario);
mkdir(plot_dir_folder)


same_low_rss = generated_samples_from_scenario(1+(3*opposite_dir_flag),experiment_name,int64(total_samples));
same_high_rss = generated_samples_from_scenario(3+(3*opposite_dir_flag),experiment_name,int64(total_samples));
same_high_rss_mat = data_cell2mat(same_high_rss);
same_low_rss_mat = data_cell2mat(same_low_rss);
[field_med_dirty,field_med_dirty_mat] = load_field_data_scenario(scenario,'med',800);
[field_low_dirty,field_low_dirty_mat] = load_field_data_scenario(scenario,'low',800);
[field_high_dirty,field_high_dirty_mat] = load_field_data_scenario(scenario,'high',800);
field_all_dirty = [field_low_dirty_mat;field_med_dirty_mat;field_high_dirty_mat];
censor_function_handle = @(x)censor_function(x,censor_noise_level,censor_pkt_size,censor_trunc_val);
[field_med_dataset_cell,field_med_per] = censor_data(field_med_dirty,censor_function_handle);
field_med_prctl = percentile_array([10,50,90],field_med_dataset_cell);
figure('Position',[1 1 2000 1000]);
total_plots = 6;
densities = linspace(0,50,total_plots);
densities = [1,5,20,22,25,27,30,40];
total_plots = length(densities);
for i = 1:total_plots
    target_density = densities(i);
    p_high = sigmoid(target_density,a_sigmoid,bias_sigmoid);
%     p_high = target_density/50;
    val_data_idx = abs(field_all_dirty(:,3)-target_density)<1;
    val_data_field = field_all_dirty(val_data_idx,1:2);
    val_data_cell = data_mat_cell(val_data_field,800);
    [val_data_cell,val_per]=censor_data(val_data_cell,censor_function_handle);
    val_data_prctl = percentile_array([10,50,90],val_data_cell);
    generated_med_cell = mix_data_cell_array(same_high_rss,same_low_rss,p_high);
    [generated_med_cell_trunc,gen_per] = censor_data(generated_med_cell,censor_function_handle);
    gen_trunc = percentile_array([10,50,90],generated_med_cell_trunc);
    subplot(2,total_plots,i);plot(1:800,gen_trunc,1:800,val_data_prctl);grid on;legend('10% gen','50% gen','90% gen','10% field','50% field','90% field');title(sprintf('d:%g p:%g',target_density,p_high));xlim([0 800]);
    subplot(2,total_plots,i+total_plots);plot(1:800,val_per,1:800,gen_per);grid on;legend(sprintf('Field Density %g PER',target_density),'Generated PER');ylabel('PER');xlabel('TxRxRange(m)');xlim([0,800]);
    
end
suptitle(sprintf('%s Mixture Interpolation a:%g b:%g',upper(scenario), a_sigmoid,bias_sigmoid));
saveas(gcf,fullfile(plot_dir_folder,sprintf('%d.jpg',total_plots)));
