function [x, P] = run_kalman_filter(x, A, P, Q, H, R, z)
    % Prediction step
    x_pred = A * x;
    P_pred = A * P * A' + Q;

    % Kalman gain calculation
    K = P_pred * H' / (H * P_pred * H' + R);

    % Update step
    x = x_pred + K * (z - H * x_pred);
    P = (eye(size(x, 1)) - K * H) * P_pred;
end