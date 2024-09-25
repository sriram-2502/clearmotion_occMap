function R = build_measurement_noise_matrix(state_size, R_occ, R_wheel)
    R = [R_occ * eye(state_size), zeros(state_size, 2);  % Noise for occupancy data
         zeros(2, state_size), R_wheel * eye(2)];  % Noise for wheel height data
end