function weights = fit_rbf(estimated_height_map, bin_params, occ_dist, centers, sigma, query_point)
   
    % Unpack binning parameters
    new_rows = bin_params.new_rows;
    new_cols = bin_params.new_cols;
    show_plots = true;
   
    %% Normalize the estimated height map
    max_val_est = max(estimated_height_map, [], "all");
    normalized_height_map = (estimated_height_map ./ max_val_est);
    
    %% Fit RBFs to the height map
    weights = get_rbf_weights(normalized_height_map, centers, sigma);
    
    %% Generate fine grid for reconstructing the height map
    [X,Y] = meshgrid(linspace(1, new_cols, new_cols), linspace(1, new_rows, new_rows));
    coordinates = [X(:), Y(:)];
    
    % Reconstruct the height map on the fine grid using the RBFs and weights
    num_centers = size(centers, 1);
    rbf_matrix = zeros(size(coordinates, 1), num_centers);
    for i = 1:num_centers
        diff = coordinates - centers(i, :);
        rbf_matrix(:, i) = exp(-sum(diff.^2, 2) / (2 * sigma^2));
    end
    
    % Reconstruct the height map on the fine grid
    fitted_height_map = reshape(rbf_matrix * weights, size(X));
    
    %% Query a height value at a specific point (optional)
    if nargin > 5
        % Query the height at the given query point (x, y)
        x_query = query_point(1);
        y_query = query_point(2);
        height_value_query = query_rbf_height([x_query, y_query], weights, centers, sigma);
    else
        height_value_query = [];  % If no query point is provided
    end
    
    %% Plot the original and fitted height maps
    if(show_plots)
        figure(11);
        
        % Plot the original estimated height map
        [x, y] = meshgrid(1:new_cols, 1:new_rows);  % Original grid
        x = x + occ_dist;
        subplot(1, 2, 1);
        surf(x, y, normalized_height_map);
        title('Original Estimated Height Map');
        xlabel('X');
        ylabel('Y');
        zlabel('Height');
        view(2)
        
        % Plot the fitted height map on the fine mesh grid
        subplot(1, 2, 2);
        surf(x, y, fitted_height_map, 'Edgecolor', 'none');
        title('Fitted Estimated Height Map on Fine Mesh Grid');
        xlabel('X');
        ylabel('Y');
        zlabel('Height');
        view(2);
        
        % Make the figure look nicer
        colormap jet;
        colorbar;
        % pause(0.2)
    end
end