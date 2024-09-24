clc; clear; close all;

%%
sID = 77859; offset = 124.143;
s = loadSession(sID, 'includeUhfd', 1);
filename = "occMsgs_"+string(sID)+".json";
fileID = fopen(filename, 'r');raw = fread(fileID, inf);
fclose(fileID);
str = char(raw');
data = jsondecode(str);

% %%
% Logs containing height curves:
% > Session 77817: timeOffsetShift = 20.058
% > Session 77818: timeOffsetShift = 111.618
% > Session 77819: timeOffsetShift = 259.489
% > Session 77821: timeOffsetShift = 312.308
% Logs containing occupancy grids:
% > Session 77824: timeOffsetShift = 52.006
% > Session 77826: timeOffsetShift = 302.973


%% 
% Compute the session distance when occ was generated
ax = gca;
for ct = 1:numel(data)
    occ = data(ct).data.occupied_upper_bound;
    data(ct).dist = interp1(s.time_offset, s.dist, data(ct).timeOffset-data(ct).timestampDelta-offset);
    occ2 = reshape(occ,80,100);
    mesh(ax, occ2)
    zlim([-0.1 0.1])
    ax.Title.String = string(data(ct).timeOffset - data(ct).timestampDelta -offset);
    disp("Progress " + string(ct/numel(data)))
    %pause(0.01)
end


%%
ct = 452 % Inspect one snapshot
data(ct).dist = interp1(s.time_offset, s.dist, data(ct).timeOffset-data(ct).timestampDelta-offset);
data(ct).dist
idx = s.dist >= data(ct).dist+5 & s.dist <= data(ct).dist+10;

figure()
subplot(212)
% hold on

occ = data(ct).data.occupied_upper_bound;
occ2 = reshape(occ,80,100);
X = 5:0.05:9.95;
Y = -2.0:0.05:1.95;
X = X + data(ct).dist;
[XX,YY] = meshgrid(X,Y);
mesh(XX, YY, occ2)

xn = 40-15;
rightRoad = occ2(xn,:);
rightRoad_ = rightRoad;
rightRoad(rightRoad==0) = NaN;
rightRoad = fillmissing(rightRoad, 'linear');
win = hann(3);
rightRoad = conv(rightRoad, win, 'same')/sum(win);
%droad = gradient(rightRoad(:))/0.05;
droad = conv(rightRoad(:), [.5 0 -.5], 'same')/0.05;
%droad = conv(droad, win, 'same')/sum(win);

subplot(211)
plot(s.dist(idx)-0*0.40, s.frBlendedRoadVel(idx)./s.speed(idx), ...
     XX(20,:), [rightRoad(:) droad], 'LineWidth', 1.5)
xlim(data(ct).dist+[5 10]+0.2)
ylim([-1,1])
grid minor
title(['Session:'+string(sID)+', Time='+string(data(ct).timeOffset-data(ct).timestampDelta-offset)+ 's'])
legend('Road Motion Slope', 'PhiGent Height (fixed)',  'PhiGent Slope')
