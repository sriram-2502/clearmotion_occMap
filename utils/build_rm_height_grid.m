function rm_height_grid = build_rm_height_grid(flHtClosest, frHtClosest, bin_params)
    % Unpack binning parameters
    M = bin_params.M;  % Number of rows in the grid
    N = bin_params.N;  % Number of columns in the grid
    wheel_offset = bin_params.wheel_offset;
    num_cells_in_bin = bin_params.num_cells_in_bin;

    % Initialize a 2D grid with zeros
    rm_height_grid = zeros(N, M);

    % Determine Y-axis (row) range for the 2D grid
    Y_range = linspace(-2, 2, N);  % Y-axis range from -2m to 2m

    % Find the row indices for the front-left and front-right wheel positions
    [~, left_wheel_idx] = min(abs(Y_range + wheel_offset));  % Front-left wheel (negative offset)
    [~, right_wheel_idx] = min(abs(Y_range - wheel_offset));  % Front-right wheel (positive offset)

    % Place the 1D height vectors (flHtClosest and frHtClosest) into the 2D grid
    tyre_patch_offset = floor(num_cells_in_bin/2);
    start_idx_left = left_wheel_idx-tyre_patch_offset; 
    end_idx_left = left_wheel_idx+tyre_patch_offset;
   
    rm_height_grid(start_idx_left:end_idx_left, :) = repmat(flHtClosest,num_cells_in_bin,1);  % Place front-left heights
    start_idx_right = right_wheel_idx-tyre_patch_offset; 
    end_idx_right = right_wheel_idx+tyre_patch_offset;
    rm_height_grid(start_idx_right:end_idx_right, :) = repmat(frHtClosest,num_cells_in_bin,1);  % Place front-right heights
end