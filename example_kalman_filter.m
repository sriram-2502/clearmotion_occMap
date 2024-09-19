% Kalman Filter for 1D position and velocity estimation with two sensors
clc; clear;

% Time step (delta t)
dt = 1;  % 1 second

% Initial State [position; velocity]
x = [0; 1];  % Initial position = 0, initial velocity = 1

% State transition matrix (A)
A = [1, dt; 
     0, 1];

% Control matrix (B) - assuming no control input
B = [0; 0];

% Measurement matrix (H) - Two sensors measuring the position
H = [1, 0;  % Sensor 1 measures position
     1, 0]; % Sensor 2 also measures position

% Process noise covariance (Q)
Q = [1, 0; 
     0, 1] * 0.001;  % Small process noise

% Measurement noise covariance (R) for two sensors
R = [2, 0;  % Sensor 1 has variance 2
     0, 5]; % Sensor 2 has variance 5

% Initial estimation error covariance (P)
P = eye(2);  % Initial uncertainty

% Number of time steps
N = 50;

% Storage for results
estimated_positions = zeros(1, N);
estimated_velocities = zeros(1, N);
true_positions = zeros(1, N);
sensor1_measurements = zeros(1, N);
sensor2_measurements = zeros(1, N);

% Simulate the system and run Kalman Filter
for k = 1:N
    % True position and velocity (simulating a real system)
    true_position = 0.5 * k^2;  % Assume true motion: position = 0.5 * t^2
    velocity = k;  % Assume velocity = t

    % Simulated noisy measurements from both sensors
    sensor1_measurement = true_position + sqrt(R(1,1)) * randn;  % Sensor 1 measurement
    sensor2_measurement = true_position + sqrt(R(2,2)) * randn;  % Sensor 2 measurement
    
    sensor1_measurements(k) = sensor1_measurement;
    sensor2_measurements(k) = sensor2_measurement;

    % Combined measurement vector
    z = [sensor1_measurement; sensor2_measurement];

    % Prediction step
    x_pred = A * x;            % Predict the next state
    P_pred = A * P * A' + Q;   % Predict the next error covariance

    % Kalman gain
    K = P_pred * H' / (H * P_pred * H' + R);

    % Update step
    x = x_pred + K * (z - H * x_pred);  % Update the state with measurements
    P = (eye(2) - K * H) * P_pred;      % Update the error covariance

    % Store the estimated position and velocity
    estimated_positions(k) = x(1);
    estimated_velocities(k) = x(2);
    true_positions(k) = true_position;
end

% Plot the results
time = 1:N;

figure;
subplot(2,1,1);
plot(time, true_positions, 'g-', 'LineWidth', 2); hold on;
plot(time, sensor1_measurements, 'rx', 'LineWidth', 1.5);
plot(time, sensor2_measurements, 'bx', 'LineWidth', 1.5);
plot(time, estimated_positions, 'b-', 'LineWidth', 2);
legend('True Position', 'Sensor 1 Measurement', 'Sensor 2 Measurement', 'Estimated Position');
xlabel('Time');
ylabel('Position');
title('Position Estimation with Two Sensors');
grid on;

subplot(2,1,2);
plot(time, estimated_velocities, 'r-', 'LineWidth', 2);
xlabel('Time');
ylabel('Velocity');
title('Estimated Velocity');
grid on;
