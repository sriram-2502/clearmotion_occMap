% Kalman Filter for 2D Grid (occupancy/height map)
clc; clear; close all;

%% load occ data

sID = 77859; offset = 124.143;
s = loadSession(sID, 'includeUhfd', 1)
filename = "occMsgs_"+string(sID)+".json";
fileID = fopen(filename, 'r');
raw = fread(fileID, inf);
fclose(fileID);
str = char(raw');
data = jsondecode(str);

% Compute the session distance when occ was generated
ax = gca;
for ct = 1:numel(data)
    occ = data(ct).data.occupied_upper_bound;
    data(ct).dist = interp1(s.time_offset, s.dist, data(ct).timeOffset-data(ct).timestampDelta-offset);
    occ2 = reshape(occ,80,100);
    mesh(ax, occ2)
    zlim([-0.1 0.1])
    ax.Title.String = string(data(ct).timeOffset - data(ct).timestampDelta -offset);
    disp("Progress " + string(ct/numel(data)))
    %pause(0.01)
end

%% test kalman filter

% Number of time steps (replace this with the actual number of steps in your data)
time_steps = 5;

% Grid dimensions (N x M) based on your example
N = 80;  % Number of rows in the grid
M = 100; % Number of columns in the grid
grid_size = N * M;

% Time step (delta t)
dt = 1;

% Initial state: height or occupancy values for the grid (all zeros to start)
x = zeros(grid_size, 1);

% State transition matrix (A)
A = eye(grid_size);  % Assuming heights/occupancies are static over time

% Measurement matrix (H)
H = eye(grid_size);  % Direct measurements for the whole grid

% Process noise covariance (Q)
Q = eye(grid_size) * 0.001;  % Small process noise (assuming little change in heights/occupancies)

% Measurement noise covariance (R)
R = eye(grid_size) * 10;  % High measurement noise (adjust based on your actual data)

% Initial estimation error covariance (P)
P = eye(grid_size);  % Initial uncertainty in the grid state

% Storage for results
estimated_heights = zeros(grid_size, time_steps);
measurements = zeros(grid_size, time_steps);

% Load and run Kalman Filter for each time step
for ct = 1:time_steps
    ct
    % Load grid data for this time step
    occ = data(ct).data.occupied_upper_bound;  % Load occupancy data
    occ2 = reshape(occ, N, M);  % Reshape to N x M grid (80x100 as per your example)
    
    % Simulated noisy measurement of the grid
    measurement = occ2(:);  % Flatten the grid into a column vector
    measurements(:, ct) = measurement;
    
    % Prediction step (Kalman filter)
    x_pred = A * x;            % Predict the next state (height or occupancy)
    P_pred = A * P * A' + Q;   % Predict the next error covariance
    
    % Kalman gain
    K = P_pred * H' / (H * P_pred * H' + R);
    
    % Update step
    x = x_pred + K * (measurement - H * x_pred);  % Update the state with measurement
    P = (eye(grid_size) - K * H) * P_pred;        % Update the error covariance
    
    % Store the estimated heights or occupancies
    estimated_heights(:, ct) = x;
end

%% Visualize the results for the last time step (you can adjust this to any other time step)
time_step_to_plot = time_steps;

% Reshape the estimated heights/occupancies back to 2D grid
estimated_heights_2d = reshape(estimated_heights(:, time_step_to_plot), N, M);
measured_heights_2d = reshape(measurements(:, time_step_to_plot), N, M);

% Plot the original measurement and estimated heights (occupancy grid)
figure;

subplot(1, 2, 1);
X = 5:0.05:9.95;
Y = -2.0:0.05:1.95;
X = X + data(ct).dist;
[XX,YY] = meshgrid(X,Y);
mesh(XX, YY, measured_heights_2d);  % Use your X, Y grid if needed
title('Measured Heights/Occupancies');
xlabel('X'); ylabel('Y'); zlabel('Height/Occupancy');
axis tight; view([140 35]);

subplot(1, 2, 2);
mesh(XX, YY, estimated_heights_2d);  % Use your X, Y grid if needed
title('Estimated Heights/Occupancies (Kalman Filter)');
xlabel('X'); ylabel('Y'); zlabel('Height/Occupancy');
axis tight; view([140 35]);