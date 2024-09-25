clc; clear; close all;
addpath utils\
addpath kalman_filter\

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
speed = s.speed;  % Vehicle speed
distances = s.dist;

% Compute the session distance when occ was generated
for ct = 1:numel(data)
    data(ct).dist = interp1(s.time_offset, s.dist, data(ct).timeOffset - data(ct).timestampDelta - offset);
end

%% Parameter setup
track_width = 1.5; % trackwidth of Q7
bin_size_in_inches = 5;  % size of a typical tyre contact patch
occ_resolution = 0.05;  % occ grid resolution (0.05 meters)

% Grid dimensions of original occupancy map
N = 80;  % Number of rows
M = 100; % Number of columns

bin_params = get_bin_params(N, M, bin_size_in_inches, track_width, occ_resolution);
binned_grid_size = bin_params.binned_grid_size;

%% Kalman Filter Parameters
start_idx = 100;
end_idx = 200;
length = end_idx - start_idx;  % Number of time steps

% Initialize state vector
state_size = binned_grid_size;  % Occupancy map size only
x = zeros(state_size, 1);  % Initial state vector: occupancy

% State transition matrix
A = eye(state_size);  % State transition matrix

% Process noise covariance (Q)
Q = eye(state_size) * 1;  % Process noise

% Initial estimation error covariance (P)
P = eye(state_size);  % Initial uncertainty

% Measurement noise covariances for the fused sources
R_occ = 10;  % Higher noise for occupancy data
R_wheel = 1;  % Lower noise for wheel height data

% Storage for results
estimated_heights = zeros(binned_grid_size, length);
measurements = zeros(state_size, length);

%% Run the Kalman Filter for each time step
for ct = start_idx:end_idx
    % Load the original occupancy data for this time step
    occ = data(ct).data.occupied_upper_bound;  % Load occupancy data
    occ2 = reshape(occ, N, M);  % Reshape to original 80x100 grid

    % Retrieve road height data for this time step
    idx = s.dist >= data(ct).dist+5 & s.dist <= data(ct).dist+10;
    front_left_ht = flRoadHt(idx)./speed(idx);  % Front-left wheel heights
    front_right_ht = frRoadHt(idx)./speed(idx);  % Front-right wheel height

    % Retrieve and bin the road heights measured at left and right wheels
    [flHtClosest, frHtClosest, rm_height_grid] = bin_road_heights(flRoadHt, frRoadHt, speed, distances, idx, bin_params);

    % Bin the occupancy map and RM heights
    occ_binned = bin_grid_map(occ2, bin_params);
    rm_height_binned = bin_grid_map(rm_height_grid, bin_params);
    
    % Combine the two sources of measurements
    combined_measurement = [occ_binned(:); flHtClosest; frHtClosest];

    % Build measurement matrix H
    H = build_measurement_matrix(state_size, left_wheel_idx, right_wheel_idx);

    % Build noise matrix R for the fused data
    R = build_measurement_noise_matrix(state_size, R_occ, R_wheel);

    % Run the Kalman filter for this time step
    [x, P] = run_kalman_filter(x, A, P, Q, H, R, combined_measurement);

    % Store the estimated heights
    estimated_heights(:, ct - start_idx + 1) = x;

    % Visualization (optional)
    visualize_heights(estimated_heights, new_rows, new_cols, data(ct).dist, ct);
end