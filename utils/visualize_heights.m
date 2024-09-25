function visualize_heights(estimated_heights, new_rows, new_cols, dist, ct)
    % Reshape the estimated heights
    estimated_heights_2d = reshape(estimated_heights(:, ct), new_rows, new_cols);

    % Create X and Y grid for visualization
    X = 5:0.05:9.95;  % Original X grid (in meters)
    Y = -2.0:0.05:1.95;  % Original Y grid (in meters)
    X = X + dist;

    % Binned X and Y grid
    new_X = linspace(min(X), max(X), new_cols);
    new_Y = linspace(min(Y), max(Y), new_rows);
    
    [XX, YY] = meshgrid(new_X, new_Y);

    % Plot the Kalman filtered estimated heights
    figure(1);
    surf(XX, YY, estimated_heights_2d);
    title('Estimated Heights/Occupancies (Kalman Filter)');
    xlabel('X (meters)');
    ylabel('Y (meters)');
    zlabel('Height/Occupancy');
    axis tight; view(2); colorbar;
end