function [flHtClosest, frHtClosest, rm_height_grid] = bin_road_heights(flRoadHt, frRoadHt, speed, distances, idx, bin_params)
    
    % Unpack bin params
    occ_resolution = bin_params.occ_resolution; % resolution of occ map
    M = bin_params.M;  % Number of columns (grid size)

    % Filter road heights based on the relevant distance
    front_left_ht = flRoadHt(idx) ./ speed(idx);  % Front-left wheel heights (normalized by speed)
    front_right_ht = frRoadHt(idx) ./ speed(idx);  % Front-right wheel height (normalized by speed)
    distance = distances(idx);  % Filtered distances

    % Create a grid of X-positions along the occupancy map (from 0 to M * resolution)
    X_grid = (0:(M-1)) * occ_resolution;

    % Interpolate road heights to match the X-grid positions
    flHtClosest = interp1(distance, front_left_ht, X_grid, 'linear', 'extrap');  % Front-left wheel height vector
    frHtClosest = interp1(distance, front_right_ht, X_grid, 'linear', 'extrap');  % Front-right wheel height vector
    rm_height_grid = build_rm_height_grid(flHtClosest, frHtClosest, bin_params);
end
