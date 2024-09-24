clc; clear; close all;

% Conversion factor
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

% Load road profile data
flRoadHt = s.flBlendedRoadVel;  % Front-left road profile 
frRoadHt = s.frBlendedRoadVel;  % Front-right road profile
speed = s.speed;
distances = s.dist;

% Compute the session distance when occ was generated
for ct = 1:numel(data)
    data(ct).dist = interp1(s.time_offset, s.dist, data(ct).timeOffset - data(ct).timestampDelta - offset);
end

%% Binning Parameters
bin_size = 5 * inch_to_meters;  % Bin size in meters
resolution = 0.05;  % Original grid resolution (0.05 meters)
num_cells_in_bin = ceil(bin_size / resolution);  % Number of cells in one bin

% Grid dimensions of original occupancy map
N = 80;  % Number of rows
M = 100; % Number of columns
mapResolution = 0.05; % 5 cm grids-

% Track width (distance between the left and right wheels)
track_width = 1.5;  % Meters

% Determine the Y positions (left-right) of the front-left and front-right wheels
wheel_offset = track_width / 2;  % Distance from the centerline to the left or right wheel

% New dimensions after binning
new_rows = floor(N / num_cells_in_bin);
new_cols = floor(M / num_cells_in_bin);
binned_grid_size = new_rows * new_cols;  % Total grid size after binning

%% Kalman Filter Parameters
start_idx = 100;
end_idx = 200;
length = end_idx - start_idx;  % Number of time steps

% Extend state vector to include 1D acceleration data
state_size = binned_grid_size;  % Occupancy map size only
x = zeros(state_size, 1);  % Initial state vector: occupancy

% State transition and measurement matrices (static for now)
A = eye(state_size);  % State transition matrix
H = eye(state_size);  % Measurement matrix

% Process noise covariance (Q) and measurement noise covariance (R)
Q = eye(state_size) * 1;  % Process noise
R = eye(state_size) * 1;  % Measurement noise

% Initial estimation error covariance (P)
P = eye(state_size);  % Initial uncertainty

% Storage for results
estimated_heights = zeros(binned_grid_size, length);
measurements = zeros(state_size, length);

%% Run the Kalman Filter for each time step
for ct = start_idx:end_idx
    % Load the original occupancy data for this time step
    occ = data(ct).data.occupied_upper_bound;  % Load occupancy data
    occ2 = reshape(occ, N, M);  % Reshape to original 80x100 grid

    % Determine where the left and right wheels are on the grid (Y-axis)
    % Assuming the vehicle centerline is in the middle of the Y-axis
    Y_range = linspace(-2, 2, new_rows);  % Y-axis range (in meters, for binned grid)
    
    % Find the nearest indices for the front-left and front-right wheel positions
    [~, left_wheel_idx] = min(abs(Y_range + wheel_offset));  % Front-left wheel (negative offset)
    [~, right_wheel_idx] = min(abs(Y_range - wheel_offset));  % Front-right wheel (positive offset)

    % Retrieve road height data for this time step
    idx = s.dist >= data(ct).dist+5 & s.dist <= data(ct).dist+10;
    front_left_ht = flRoadHt(idx)./speed(idx);  % Front-left wheel heights
    front_right_ht = frRoadHt(idx)./speed(idx);  % Front-right wheel height

    % match road heights with occ distances
    for i = 1:M
        currentX = (i - 1) * mapResolution;
        [~, closestIdx] = min(abs(distances - currentX));
        flHtClosest = front_left_ht(closestIdx);
        frHtClosest = front_right_ht(closestIdx);
    end

    % Add the road height data to the corresponding grid cells
    % Left wheel (modify the binned occupancy grid at left_wheel_idx)
    occ2(left_wheel_idx, :) = occ2(left_wheel_idx, :) + flHtClosest;

    % Right wheel (modify the binned occupancy grid at right_wheel_idx)
    occ2(right_wheel_idx, :) = occ2(right_wheel_idx, :) + frHtClosest;

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
            
            % Filter out zeros and compute the mean
            non_zero_values = sub_grid(sub_grid ~= 0);
            if ~isempty(non_zero_values)
                occ_binned(i, j) = mean(non_zero_values);
            else
                occ_binned(i, j) = 0;  % Set to zero if no non-zero values
            end
        end
    end

    % Flatten the binned occupancy grid into a vector for Kalman filter
    occ_binned_flat = occ_binned(:);

    % Store the noisy measurement (occupancy grid with wheel accelerations added)
    measurements(:, ct - start_idx + 1) = occ_binned_flat;

    % Kalman filter prediction step
    x_pred = A * x;  % Predict the next state
    P_pred = A * P * A' + Q;  % Predict the next error covariance

    % Kalman gain calculation
    K = P_pred * H' / (H * P_pred * H' + R);

    % Kalman filter update step
    x = x_pred + K * (occ_binned_flat - H * x_pred);  % Update the state with the new measurement
    P = (eye(state_size) - K * H) * P_pred;  % Update the error covariance

    % Store the estimated heights/occupancies
    estimated_heights(:, ct - start_idx + 1) = x;

    % Reshape the estimated heights/occupancies back to 2D grid
    estimated_heights_2d = reshape(estimated_heights(:, ct - start_idx + 1), new_rows, new_cols);
    measured_heights_2d = reshape(measurements(:, ct - start_idx + 1), new_rows, new_cols);
    
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
