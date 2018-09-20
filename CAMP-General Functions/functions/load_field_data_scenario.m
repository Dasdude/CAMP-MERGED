function [dataset_cell_dirty,dataset_mat_dirty,density_east,density_west,density_same] = load_field_data_scenario(scenario,density_level,d_max,per_included)
%LOAD_FIELD_DATA_SCENARIO Summary of this function goes here
%  'same/opposite','low/med/high'
same_low_up = 15;
same_med_up = 30;
opposite_low_up = 15;
opposite_med_up = 30;
file_string = sprintf('%s/s %d %d o %d %d/data results/%s %s.csv','Seperated DensityPER',same_low_up,same_med_up,opposite_low_up,opposite_med_up,scenario,density_level);
input  = file_string;
% file_name_string = sprintf('%s/%s Direction %s Density %d to %d',experiment_name,mode{1},mode{2},mode{3},mode{4});
csv_data = readtable(input,'ReadVariableNames',true);
if ~per_included
    csv_data = csv_data(csv_data.RSS>-100,:);
end
dataset_mat_dirty = [csv_data.Range,csv_data.RSS,csv_data.average_density];
dataset_cell_dirty = data_mat_cell(dataset_mat_dirty,d_max);

if strcmp(scenario,'same')
    density_west = table2array(csv_data(:,28:35));
    density_east = table2array(csv_data(:,18:25));
    density = (density_east.*repmat(1-csv_data.ego_side,1,8))+(density_west.*repmat(csv_data.ego_side,1,8));
    density_same = array2table(density);
else
    density_same = [];
end
density_west = csv_data(:,28:35);
density_east = csv_data(:,18:25);
% dataset_cell_dirty = data_mat_cell(dataset_mat_dirty,d_max);
end

