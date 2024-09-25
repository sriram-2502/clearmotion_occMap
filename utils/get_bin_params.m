function bin_params = get_bin_params(N, M, bin_size_in_inches, track_width, occ_resolution)
    % Conversion factor
    inch_to_meters = 0.0254;

    % Convert bin size from inches to meters
    bin_size = bin_size_in_inches * inch_to_meters;

    % Number of cells in one bin
    num_cells_in_bin = ceil(bin_size / occ_resolution);

    % Calculate the Y-offset (distance between the left and right wheels)
    wheel_offset = track_width / 2;

    % Calculate the new grid dimensions after binning
    new_rows = floor(N / num_cells_in_bin);
    new_cols = floor(M / num_cells_in_bin);
    binned_grid_size = new_rows * new_cols;

    % Populate the struct with all the calculated parameters
    bin_params = struct();
    bin_params.bin_size = bin_size;
    bin_params.num_cells_in_bin = num_cells_in_bin;
    bin_params.N = N;  % Number of rows in the original grid
    bin_params.M = M;  % Number of columns in the original grid
    bin_params.occ_resolution = occ_resolution;
    bin_params.track_width = track_width;
    bin_params.wheel_offset = wheel_offset;
    bin_params.new_rows = new_rows;
    bin_params.new_cols = new_cols;
    bin_params.binned_grid_size = binned_grid_size;
end