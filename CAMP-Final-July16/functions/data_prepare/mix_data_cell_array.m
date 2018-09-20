function [mixed_datacell] = mix_data_cell_array(array_x,array_y,ratio_x)
%MIX_DATA_CELL_ARRAY Summary of this function goes here
%   Detailed explanation goes here
mixed_datacell = cell(1,length(array_x));
    for i = 1:length(array_x)
        total_samples = length(array_x{i});
        x_total = int64(total_samples*ratio_x);
        y_total = int64(total_samples*(1-ratio_x));
        x_idx = randsample(total_samples,x_total);
        y_idx = randsample(total_samples,y_total);
        x_sel_array = array_x{i}(x_idx);
        y_sel_array = array_y{i}(y_idx);
        mixed_datacell{i} = [x_sel_array;y_sel_array];
    end
end

