clc; clear; close all;
inch_to_meters = 0.0254;

%% Load occupancy data
sID = 77859; offset = 124.143;
s = loadSession(sID, 'includeUhfd', 1);
filename = "occMsgs_" + string(sID) + ".json";
fileID = fopen(filename, 'r');
raw = fread(fileID, inf);
fclose(fileID);
str = char(raw');
data = jsondecode(str);

% Compute the session distance when occ was generated
for ct = 1:numel(data)
    data(ct).dist = interp1(s.time_offset, s.dist, data(ct).timeOffset - data(ct).timestampDelta - offset);
end

%% Binning Parameters
bin_size = 5*inch_to_meters;  % Bin size in meters
resolution = 0.05;  % Original grid resolution (0.05 meters)
num_cells_in_bin = ceil(bin_size / resolution);  % Number of cells in one bin

% Grid dimensions of original occupancy map
N = 80;  % Number of rows
M = 100; % Number of columns

% New dimensions after binning
new_rows = floor(N / num_cells_in_bin);
new_cols = floor(M / num_cells_in_bin);
binned_grid_size = new_rows * new_cols;  % Total grid size after binning

%% Kalman Filter Parameters
% set the start and end points
start_idx = 100;
end_idx = 200;
length = end_idx-start_idx;  % Number windows to store

% Initial state for binned data (all zeros)
x = zeros(binned_grid_size, 1);

% State transition and measurement matrices (assume static for both)
A = eye(binned_grid_size);  % State transition matrix (static)
H = eye(binned_grid_size);  % Measurement matrix (static)

% Process noise covariance (Q) and measurement noise covariance (R)
Q = eye(binned_grid_size) * 1;  % Small process noise
R = eye(binned_grid_size) * 1;     % Larger measurement noise

% Initial estimation error covariance (P)
P = eye(binned_grid_size);  % Initial uncertainty

% Storage for results
estimated_heights = zeros(binned_grid_size, length);
measurements = zeros(binned_grid_size, length);

%% Run the Kalman Filter for each time step
for ct = start_idx:end_idx
    % Load the original occupancy data for this time step
    occ = data(ct).data.occupied_upper_bound;  % Load occupancy data
    occ2 = reshape(occ, N, M);  % Reshape to original 80x100 grid

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

    % Flatten the binned grid into a vector (for Kalman filter processing)
    measurement = occ_binned(:);
    measurements(:, ct) = measurement;  % Store the noisy measurement

    % Kalman filter prediction step
    x_pred = A * x;  % Predict the next state
    P_pred = A * P * A' + Q;  % Predict the next error covariance

    % Kalman gain calculation
    K = P_pred * H' / (H * P_pred * H' + R);

    % Kalman filter update step
    x = x_pred + K * (measurement - H * x_pred);  % Update the state with the new measurement
    P = (eye(binned_grid_size) - K * H) * P_pred;  % Update the error covariance

    % Store the estimated heights/occupancies
    estimated_heights(:, ct) = x;
   
    % Reshape the estimated heights/occupancies back to 2D grid
    estimated_heights_2d = reshape(estimated_heights(:, ct), new_rows, new_cols);
    measured_heights_2d = reshape(measurements(:, ct), new_rows, new_cols);
    
    % Create X and Y grid for visualization (based on binned dimensions)
    X = 5:0.05:9.95;  % Original X grid (in meters)
    Y = -2.0:0.05:1.95;  % Original Y grid (in meters)
    X = X + data(ct).dist;
    
    % Binned X and Y grid (resized to match the binned grid dimensions)
    new_X = linspace(min(X), max(X), new_cols);
    new_Y = linspace(min(Y), max(Y), new_rows);
    
    % Plot both the original measured data (binned) and the Kalman filtered data
    figure(1);
    
    % Plot the binned measured data
    subplot(2, 1, 1);
    [XX, YY] = meshgrid(new_X, new_Y);
    surf(XX, YY, measured_heights_2d);  % Binned measured heights
    title('Binned Measured Heights/Occupancies');
    xlabel('X (meters)');
    ylabel('Y (meters)');
    zlabel('Height/Occupancy');
    axis tight; view(2); colorbar;
    
    % Plot the estimated data from Kalman filter
    subplot(2, 1, 2);
    surf(XX, YY, estimated_heights_2d);  % Kalman-filtered estimated heights
    title('Estimated Heights/Occupancies (Kalman Filter)');
    xlabel('X (meters)');
    ylabel('Y (meters)');
    zlabel('Height/Occupancy');
    axis tight; view(2); colorbar;

    pause(0.5)
end