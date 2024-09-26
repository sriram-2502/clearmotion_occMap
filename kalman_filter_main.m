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
new_rows = bin_params.new_rows;
new_cols = bin_params.new_cols;

binned_grid_size = bin_params.binned_grid_size;

% storage
start_idx = 100;
end_idx = 200;
length = end_idx - start_idx;  % Number of time steps
estimated_heights = zeros(binned_grid_size, length);
measurements = zeros(binned_grid_size, length);

% rbf params
sigma = 1; grid_spacing = 1;
[center_x, center_y] = meshgrid(1:grid_spacing:new_cols, 1:grid_spacing:new_rows);
centers = [center_x(:), center_y(:)];  % Flatten the grid of centers

%% Kalman Filter Parameters

% Initialize state vector
state_size = 2*binned_grid_size;  % occ_map size + rm_heigt_map size
x = zeros(state_size, 1);  % Initial state vector: occupancy

% State transition matrix
A = eye(state_size);  % State transition matrix

% Process noise covariance (Q)
% Higher process noise for occupancy data, lower for road height data
Q_occ = 1 * eye(state_size/2);  % Process noise for occupancy
Q_rm = 0.1 * eye(state_size/2);  % Lower process noise for road heights
Q = blkdiag(Q_occ, Q_rm);  % Combine into a block diagonal matrix

% Initial estimation error covariance (P)
% Higher initial uncertainty for occupancy, lower for road height
P_occ = 1 * eye(state_size/2);  % Initial uncertainty for occupancy
P_rm = 0.1 * eye(state_size/2);  % Lower initial uncertainty for road height
P = blkdiag(P_occ, P_rm);  % Combine into a block diagonal matrix

% Measurement noise covariances for the fused sources
noise_occ = 10;  % Higher noise for occupancy data
noise_rm = 1;  % Lower noise for wheel height data

% Build the measurement matrix H
H_occ = eye(state_size/2);
H_rm = eye(state_size/2);
H = blkdiag(H_occ, H_rm);  % Block diagonal matrix to combine occupancy and road height matrices

% Build the measurement noise covariance matrix R
R_occ = noise_occ * eye(state_size/2);  % Higher noise for occupancy data
R_rm = noise_rm * eye(state_size/2);     % Lower noise for road height data
R = blkdiag(R_occ, R_rm);  % Block diagonal matrix for the combined noise

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
    combined_measurement = [occ_binned(:); rm_height_binned(:)];

    % Run the Kalman filter for this time step
    [x, P] = run_kalman_filter(x, A, P, Q, H, R, combined_measurement);

    % Store the estimated heights
    estimated_height_map = reshape(x(1:state_size/2),new_rows, new_cols);
    estimated_heights(:, ct - start_idx + 1) = x(1:state_size/2); % get only first N states

    % fit rbfs and form density function
    weights_est = fit_rbf(estimated_height_map, bin_params, data(ct).dist, centers, sigma);

    % Visualization (optional)
    % visualize_heights(estimated_height_map, data(ct).dist, bin_params);
end