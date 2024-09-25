function grid_map_binned = bin_grid_map(grid_map, bin_params)
    
    % unpack bin params
    num_cells_in_bin = bin_params.num_cells_in_bin;
    N = bin_params.N;  % Number of rows in the original grid
    M = bin_params.M;  % Number of columns in the original grid
    new_rows = bin_params.new_rows;
    new_cols = bin_params.new_cols;

    grid_map_binned = zeros(new_rows, new_cols);  % Initialize the binned occupancy map
    for i = 1:new_rows
        for j = 1:new_cols
            % Extract sub-grid
            row_start = (i-1) * num_cells_in_bin + 1;
            row_end = min(i * num_cells_in_bin, N);  % Avoid out-of-bound errors
            col_start = (j-1) * num_cells_in_bin + 1;
            col_end = min(j * num_cells_in_bin, M);
            sub_grid = grid_map(row_start:row_end, col_start:col_end);
            
            % Compute the mean of the non-zero values
            non_zero_values = sub_grid(sub_grid ~= 0);
            if ~isempty(non_zero_values)
                grid_map_binned(i, j) = mean(non_zero_values);
            end
        end
    end
end