clc; clear; close all;
save_height_map = true;
save_accel_map = true;


%% load session and set params
sID = 77859; offset = 124.143;

% Load session data
s = loadSession(sID, 'includeUhfd', 1);

% Load occupancy message data
filename = "occMsgs_"+string(sID)+".json";
fileID = fopen(filename, 'r');
raw = fread(fileID, inf);
fclose(fileID);
str = char(raw');
data = jsondecode(str);

% Parameters
mapResolution = 0.05;    % Resolution of the map (meters per cell)
mapSizeY = 80;           % Fixed Y size (same as occupancy map size)
dt = 0.02;               % update based on session data
trackWidth = 1.5;        % Track width in meters
gridSizeX = 100;         % Number of grid columns
gridSizeY = 80;          % Number of grid rows

% Set the global origin for Y axis
originY = mapSizeY / 2;

% find trackwidth locations
% Find the grid rows for the left and right wheel positions
Y = -2.0:mapResolution:1.95; % Y limits
[~, leftWheelRow] = min(abs(Y - (-trackWidth / 2)));
[~, rightWheelRow] = min(abs(Y - (trackWidth / 2)));

% set the start and end points
start_idx = 800;
end_idx = 900;

%% For initializatinon, set static map and static accel grid as the current data
% Calculate the distance corresponding to this occupancy grid
ct = start_idx;
data(ct).dist = interp1(s.time_offset, s.dist, data(ct).timeOffset - data(ct).timestampDelta - offset);

% Extract occupancy grid and reshape
occ = data(ct).data.occupied_upper_bound;
occ2 = reshape(occ, gridSizeY, gridSizeX);  % Use same grid size as static map
heightMap = occ2;

% Initialize acceleration grid for this specific snapshot
accelMap = nan(size(occ2));
distances = s.dist;
flAccelRaw = s.flWheelAccel;
frAccelRaw = s.frWheelAccel;

% Fill acceleration grid based on current distances
for i = 1:gridSizeX
    currentX = (i - 1) * mapResolution;
    [~, closestIdx] = min(abs(distances - currentX));
    flAccelClosest = flAccelRaw(closestIdx);
    frAccelClosest = frAccelRaw(closestIdx);
    accelMap(leftWheelRow, i) = flAccelClosest;
    accelMap(rightWheelRow, i) = frAccelClosest;

    % set zero accel values at the edges
    accelMap(1,i) = 0;
    accelMap(end,i) = 0;
end

% fill nan values and apply smoothing
accelMapFilled = fillmissing(accelMap, 'linear');
% smoothWindowSize = 5; % Adjust as needed
% h = fspecial('gaussian', smoothWindowSize, smoothWindowSize/6); % Gaussian filter
% accelMapSmoothed = imfilter(accelMapFilled, h, 'replicate');
accelMapSmoothed = accelMapFilled;

% % Display the results
% figure;
% subplot(2,1,1);
% imagesc(accelMap); % Display original grid with known values
% colorbar;
% title('Original Acceleration Grid with Known Values');
% 
% subplot(2,1,2);
% imagesc(accelMapSmoothed); % Display interpolated and smoothed grid
% colorbar;
% title('Interpolated and Smoothed Acceleration Grid');



%% Loop through all occupancy data and accumulate into static map and static accel grid
% Initialize the last distance tracker
lastDist = data(ct).dist;

for ct = start_idx+1:end_idx
    % Calculate the distance corresponding to this occupancy grid
    data(ct).dist = interp1(s.time_offset, s.dist, data(ct).timeOffset - data(ct).timestampDelta - offset);
    
    % Extract occupancy grid and reshape
    occ = data(ct).data.occupied_upper_bound;
    occ2 = reshape(occ, gridSizeY, gridSizeX);  % Use same grid size as static map
    occSizeX = size(occ2, 2);  % Get the num cols (X size) of the occupancy grid
    occSizeY = size(occ2, 1);  % Get the num rows (Y size) of the occupancy grid

    % Compute the global position of this grid
    globalDist = data(ct).dist;  % Distance of this occupancy grid from the start
    
    % Check if static map needs to be updated
    if globalDist > lastDist
        % Calculate the number of new columns to add
        numNewColumns = floor((globalDist - lastDist) / mapResolution);

        if numNewColumns > 0
            % Expand the static map if necessary
            if isempty(heightMap)
                heightMap = nan(gridSizeY, numNewColumns);  % Initialize with the number of new columns
            else
                % Expand static map in X direction if necessary
                heightMap = [heightMap, nan(gridSizeY, numNewColumns)];
            end
            
            % Expand static acceleration grid
            if isempty(accelMap)
                accelMap = nan(gridSizeY, numNewColumns);
            else
                accelMap = [accelMap, nan(gridSizeY, numNewColumns)];
            end
            
            % Update the static map and acceleration grid with new data
            validStartIndex = size(heightMap, 2) - numNewColumns + 1;
            validEndIndex = size(heightMap, 2);
            
            % Map the occupancy grid to the correct indices
            occUpdateLength = validEndIndex - validStartIndex + 1;
            occStartIndex = occSizeX - occUpdateLength + 1;
            occEndIndex = occSizeX;
            
            % Ensure indices are within bounds
            if validStartIndex <= size(heightMap, 2)
                heightMap(:, validStartIndex:validEndIndex) = occ2(:, occStartIndex:occEndIndex);
            end
            
            % Get the wheel accelerations for the new grid columns
            distances = s.dist; % All recorded distances
            flAccelRaw = s.flWheelAccel;  % Front-left wheel acceleration
            frAccelRaw = s.frWheelAccel;  % Front-right wheel acceleration

            for i = 1:numNewColumns
                % X position for this grid column
                currentX = (size(heightMap, 2) - numNewColumns + i) * mapResolution;

                % Find the closest distance in recorded data
                [~, closestIdx] = min(abs(distances - currentX));

                % Get the corresponding acceleration values
                flAccelClosest = flAccelRaw(closestIdx);
                frAccelClosest = frAccelRaw(closestIdx);

                % Place the closest acceleration values into the grid
                accelMap(leftWheelRow, i + validStartIndex - 1) = flAccelClosest;
                accelMap(rightWheelRow, i + validStartIndex - 1) = frAccelClosest;

                % set zero accel values at the edges
                accelMap(1, i + validStartIndex - 1) = 0;
                accelMap(end, i + validStartIndex - 1) = 0;
            end
            
            % fill nan values and apply smoothing
            accelMapFilled = fillmissing(accelMap, 'linear');
            % smoothWindowSize = 5; % Adjust as needed
            % h = fspecial('gaussian', smoothWindowSize, smoothWindowSize/6); % Gaussian filter
            % accelMapSmoothed = imfilter(accelGridFilled, h, 'replicate');

            accelMapSmoothed = accelMapFilled;
            % Update the last distance
            lastDist = globalDist;
        end
    end
   
    % Create 2D plots for the static map and local occupancy grid
    figure(1)
    
    % set colormap limits
    cmin = -0.25;
    cmax = 0.25;

    % Plot local occupancy grid
    subplot(3,4,[2,3])
    % Show local occ
    X = 5:0.05:9.95;
    Y = -2.0:0.05:1.95;
    X = X + data(ct).dist;  % Update X based on current distance
    [XX, YY] = meshgrid(X, Y);
    
    % Display local occupancy grid as a 2D image
    imagesc(X(1,:), Y(:,1), occ2); 
    colorbar; 
    clim([cmin cmax]);  

    xlabel('X [m]');
    ylabel('Y [m]');
    title('Local Occupancy Grid');
    
    % Plot global static map
    subplot(3,4,[5,6,7,8])
    
    % Show global occ
    [XStatic, YStatic] = meshgrid((0:size(heightMap, 2) - 1) * mapResolution, ...
                                  (-originY:(mapSizeY - originY - 1)) * mapResolution);
    XStatic = XStatic + data(start_idx).dist;  % Update X based on current distance

    % Display global static map as a 2D image
    imagesc(XStatic(1,:), YStatic(:,1), heightMap);  % Transpose staticMap to match XStatic and YStatic dimensions
    colorbar; 
    clim([cmin cmax]);

    xlabel('X [m]');
    ylabel('Y [m]');
    title(['Stitched Static Map - Distance: ' num2str(globalDist)]);

    % Update the plot with the filtered static map
    subplot(3,4,[9,10,11,12])
    filterSize = 10;
    
    % Apply the moving average filter
    heightMapFiltered = filter2(ones(filterSize) / filterSize^2, heightMap, 'same');

    imagesc(XStatic(1,:), YStatic(:,1), heightMapFiltered);
    clim([cmin cmax]);
    colorbar;
    xlabel('X [m]');
    ylabel('Y [m]');
    title(['Filtered Static Map - Distance: ' num2str(globalDist)]);
       
    % Draw the plots
    drawnow;
    
    % Display progress
    disp("Progress: " + string(ct) + "/" + string(numel(data)));
    
    % visualize accel grid
    figure(2);
    subplot(3,4,[2,3])
    plot(s.dist, [s.flWheelAccel, s.frWheelAccel], 'LineWidth', 1.5);
    xlim([X(1), X(end)]);
    xlabel('Distance [m]');
    ylabel('Acceleration [m/s^2]');
    title('Raw Wheel Acceleration Data');
    legend('FL wheel accel', 'FR wheel accel');

    subplot(3,4,[5,6,7,8])
    imagesc(XStatic(1,:), YStatic(:,1), accelMap);  % Display the 2D grid for the selected session
    colorbar;
    xlabel('X [m]');
    ylabel('Y [m]');
    title(['Acceleration Grid - Distance: ' num2str(globalDist)]);
    
    subplot(3,4,[9,10,11,12])
    imagesc(XStatic(1,:), YStatic(:,1), accelMapSmoothed);  % Display the 2D grid for the selected session
    colorbar;
    xlabel('X [m]');
    ylabel('Y [m]');
    title(['Interpolated and Smoothed Acceleration Grid - Distance: ' num2str(globalDist)]);

end

%% save data
if(save_height_map)
    % Save the filtered static map to a .mat file
    save('heightMapFiltered.mat', 'heightMapFiltered');
end
if(save_accel_map)
    % Save the filtered static map to a .mat file
    save('accelMapSmoothed.mat', 'accelMapSmoothed');
end