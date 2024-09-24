% Bin the occ map with avg radius of 5 inch
clc; clear; close all;
inch_to_meters = 0.0254;

%% load occ data

sID = 77859; offset = 124.143;
s = loadSession(sID, 'includeUhfd', 1)
filename = "occMsgs_"+string(sID)+".json";
fileID = fopen(filename, 'r');
raw = fread(fileID, inf);
fclose(fileID);
str = char(raw');
data = jsondecode(str);

%% Original occupancy map (occ2) dimensions: 80x100

% set the start and end points
start_idx = 100;
end_idx = 900;

% Grid dimensions of original occupancy map
N = 80;  % Number of rows
M = 100; % Number of columns

for ct = start_idx+1:end_idx
    % Calculate the distance corresponding to this occupancy grid
    data(ct).dist = interp1(s.time_offset, s.dist, data(ct).timeOffset - data(ct).timestampDelta - offset);

    % X and Y dimensions are based on your original grid definition.
    X = 5:0.05:9.95;   % X-coordinates of the grid (in meters)
    Y = -2.0:0.05:1.95; % Y-coordinates of the grid (in meters)
    X = X + data(ct).dist;
    
    % Reshape occ map into 2D grid (assuming it’s already loaded)
    occ2 = reshape(data(ct).data.occupied_upper_bound, 80, 100);
    
    % Define the bin size (5 inches) and the resolution (0.05 inches)
    bin_size = 5*inch_to_meters;  % 5 inches
    resolution = 0.05;  % Each grid cell represents 0.05 meters
    
    % Number of cells in a 5-inch bin
    num_cells_in_bin = ceil(bin_size / resolution);  % 5 inch / 0.05 meters
    
    % Calculate the new size for the binned occupancy map
    new_rows = floor(size(occ2, 1) / num_cells_in_bin);
    new_cols = floor(size(occ2, 2) / num_cells_in_bin);
    
    % Initialize binned occupancy map (occ_binned)
    occ_binned = zeros(new_rows, new_cols);
    
    % Apply binning (reduce resolution)
    occ_binned = zeros(new_rows, new_cols);  % Initialize the binned occupancy map
    for i = 1:new_rows
        for j = 1:new_cols
            % Extract the sub-grid to average
            row_start = (i-1) * num_cells_in_bin + 1;
            row_end = min(i * num_cells_in_bin, N);  % Avoid out-of-bound errors
            col_start = (j-1) * num_cells_in_bin + 1;
            col_end = min(j * num_cells_in_bin, M);
    
            % Extract sub-grid
            sub_grid = occ2(row_start:row_end, col_start:col_end);
            
            % Filter out zeros
            non_zero_values = sub_grid(sub_grid ~= 0);  % Keep only non-zero values
    
            % Compute the mean if there are non-zero values, else set to zero
            if ~isempty(non_zero_values)
                occ_binned(i, j) = mean(non_zero_values);
            else
                occ_binned(i, j) = 0;  % Set to zero if no non-zero values
            end
        end
    end
    
    % Create new X and Y for the binned occupancy map
    new_X = linspace(min(X), max(X), new_cols);
    new_Y = linspace(min(Y), max(Y), new_rows);
    
    % Plot both the original and binned occupancy maps
    figure(1);
    
    % Plot the original high-resolution occupancy map
    subplot(2, 1, 1);  % Create subplot for the original map
    surf(X, Y, occ2);  % Use the original X, Y, and occ2 (high-res)
    title('Original Occupancy Map');
    ylim([-2,2]);
    xlabel('X (meters)');
    ylabel('Y (meters)');
    zlabel('Occupancy');
    axis tight;
    view(2);  % Adjust the viewing angle for a better 3D perspective
    
    % Plot the binned (low-resolution) occupancy map
    subplot(2, 1, 2);  % Create subplot for the binned map
    surf(new_X, new_Y, occ_binned);  % Plot the binned map with reduced resolution
    title('Binned Occupancy Map');
    ylim([-2,2]);
    xlabel('X Bins (meters)');
    ylabel('Y Bins (meters)');
    zlabel('Occupancy');
    axis tight;
    view(2);  % Match the viewing angle for consistency
    pause(0.2)
end