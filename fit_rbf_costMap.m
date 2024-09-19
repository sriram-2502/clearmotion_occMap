%% fit rbfs for height maps
% Fit Gaussian RBFs to the cost map data and plot the results using a fine mesh grid.

clc; clear; close all;
resolution = 0.05;

%% load the height map and accel map

loadHeightMap = load('heightMapFiltered.mat');
heightMap_ = loadHeightMap.heightMapFiltered;

loadAccelMap = load('accelMapSmoothed.mat');
accelMap_ = loadAccelMap.accelMapSmoothed;

% find slope based on heightMap 
[gradientX, gradientY] = gradient(heightMap_);
slopeMagnitude_ = sqrt(gradientX.^2 + gradientY.^2);

% set zero to eges
paddingSizeX = 10;
paddingSizeY = 20;
% Get the size of the slopeMagnitude matrix
[numRows, numCols] = size(slopeMagnitude_);
% Set top and bottom padding to zero
slopeMagnitude_(1:paddingSizeY, :) = 0;                 % Top edge
slopeMagnitude_(numRows-paddingSizeY+1:numRows, :) = 0; % Bottom edge
% Set left and right padding to zero
slopeMagnitude_(:, 1:paddingSizeX) = 0;                 % Left edge
slopeMagnitude_(:, numCols-paddingSizeX+1:numCols) = 0; % Right edge

% normalize and square the data
max_val_ht = max(heightMap_,[],"all");
heightMap = heightMap_./max_val_ht;
heightMap = heightMap.^2;

max_val_ac = max(accelMap_,[],"all");
accelMap = accelMap_./max_val_ac;
accelMap = accelMap.^2;

maxSlope = max(slopeMagnitude_, [], 'all');
slopeMagnitude = slopeMagnitude_ ./ maxSlope;

% set relatvie weights for heights and costs
heightWeight = 1;
accelWeight = 10;
slopeWeight = 5;

costMap = heightWeight.*heightMap + accelWeight.*accelMap + slopeWeight.*slopeMagnitude;

% Get the size of the cost map
[rows, cols] = size(costMap);

% Generate the grid of coordinates
[x, y] = meshgrid(1:cols, 1:rows);
coordinates = [x(:), y(:)];

figure(1)
subplot(2,3,1)
surf(x.*resolution, y.*resolution, heightMap, 'EdgeColor','none');
subplot(2,3,2)
surf(x.*resolution, y.*resolution, slopeMagnitude, 'EdgeColor','none');
subplot(2,3,3)
surf(x.*resolution, y.*resolution, accelMap, 'EdgeColor','none');

figure(2);
% Plot the original height map
subplot(1, 2, 1);
surf(x.*resolution, y.*resolution, costMap, 'EdgeColor','none');
title('True Cost Map');
xlabel('X');
ylabel('Y');
zlabel('Cost');


%% Fit RBFs for the cost map

grid_spacing = 10;  % Adjust the spacing between RBF centers
sigma = 10;  % Tune sigma as needed

% Generate the uniform grid centers
[center_x, center_y] = meshgrid(1:grid_spacing:cols, 1:grid_spacing:rows);
centers = [center_x(:), center_y(:)];

% Fit RBFs and get the weights
weights = fit_rbf(costMap, centers, sigma);

% Generate a fine grid of coordinates
[x_fine, y_fine] = meshgrid(linspace(1, cols, cols*1), linspace(1, rows, rows*1));
coordinates_fine = [x_fine(:), y_fine(:)];

% Reconstruct the height map using the RBFs and weights on the fine grid
num_centers = size(centers, 1);
rbf_matrix_fine = zeros(size(coordinates_fine, 1), num_centers);
for i = 1:num_centers
    diff = coordinates_fine - centers(i, :);
    rbf_matrix_fine(:, i) = exp(-sum(diff.^2, 2) / (2 * sigma^2));
end
fitted_cost_map_fine = reshape(rbf_matrix_fine * weights, size(x_fine));

% Query a height value at point (x, y)
x_query = 11;  % Replace with your x-coordinate of the query point
y_query = 11;  % Replace with your y-coordinate of the query point
costValue= query_rbf_height([x_query, y_query], weights, centers, sigma);
% disp(['Height at (', num2str(x_query), ', ', num2str(y_query), ') is ', num2str(costValue)]);

% Plot the fitted cost map on a fine mesh grid
figure(2);
subplot(1, 2, 2);
surf(x_fine.*resolution, y_fine.*resolution, fitted_cost_map_fine, 'Edgecolor','none');
title('Fitted Cost Map on Fine Mesh Grid');
xlabel('X');
ylabel('Y');
zlabel('Cost');

% Make the figure look nicer
colormap jet;
colorbar;