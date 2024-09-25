function H = build_measurement_matrix(state_size)
    H = [eye(state_size);  % Identity matrix for occupancy data
         eye(state_size)];  % Identity matrix for RM height grid data
end