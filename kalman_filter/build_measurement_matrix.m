function H = build_measurement_matrix(state_size, left_wheel_idx, right_wheel_idx)
    H = [eye(state_size);  % Identity matrix for occupancy data
         zeros(2, state_size)];  % Initialize with zeros for wheel height data
    H(state_size + 1, left_wheel_idx) = 1;  % Set front-left wheel height position
    H(state_size + 2, right_wheel_idx) = 1;  % Set front-right wheel height position
end