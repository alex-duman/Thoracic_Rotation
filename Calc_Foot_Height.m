function [LheelH, RheelH, LtoeH, RtoeH, L_legL, R_legL] = Calc_Foot_Height(K,speed,Trd_time,plt,fileName)
% Function estimates the left and right heel height during standing which
% can be used to determine touchdown and liftoff times (see
% Calc_Stride_Metrics.m). A secondary function is to estimate Leg length
% for both the left and right legs to satisfy the corresponding inputs for
% Calc_Stride_Metrics.m. This function assumes that the treadmill
% accelerates at 0.4 m/s^2.

% INPUTS
% K - table containing the kinematic variables of interest
% speed - calculated speed of treadmill
% Trd_time - vector of time in reference frame of treadmill (1-2kHz)
% plt - string set to 'Y...' or 'y...' if you'd like to see a visualization
%       of the foot marker position and the period used for estimating heel 
%       height
% fileName - string for the full name of the trial

% OUTPUTS
% LheelH - left heel height (averaged over 5s period during no motion)
% RheelH - right heel height (averaged over 5s period during no motion)
% LtoeH - left toe height (averaged over 5s period during no motion)
% RtoeH - right toe height (averaged over 5s period during no motion)
% L_legL - scalar of the left leg length in millimeters
% R_legL - scalar of the right leg length in millimeters

idx_motion = find(speed >= 0.1,1,'first'); % find index where speed first reaches 0.1 m/s (1 second after treadmill starts)
t_start = Trd_time(idx_motion) - 7.5; % measuring start of standing as ~7.5s prior to start of treadmill
t_end = t_start + 5; % using 5s period to estimate heel height
idx_start = find(K.time == round(t_start,2)); % round to 2 decimals since vicon measured at 100Hz
idx_end = find(K.time == round(t_end,2));

% Calculating average heel height for both feet
LheelH = mean(K.Lhee_z(idx_start:idx_end));
RheelH = mean(K.Rhee_z(idx_start:idx_end));
LtoeH = mean(K.Ltoe_z(idx_start:idx_end));
RtoeH = mean(K.Rtoe_z(idx_start:idx_end));

% Defining Points for Calculating Leg Length
v_Lleg = [K.Lasi_x(idx_start:idx_end), K.Lasi_y(idx_start:idx_end), K.Lasi_z(idx_start:idx_end)] - [K.Lank_x(idx_start:idx_end), K.Lank_y(idx_start:idx_end), K.Lank_z(idx_start:idx_end)]; % vector from left ankle -> left ASI hip marker
v_Rleg = [K.Rasi_x(idx_start:idx_end), K.Rasi_y(idx_start:idx_end), K.Rasi_z(idx_start:idx_end)] - [K.Rank_x(idx_start:idx_end), K.Rank_y(idx_start:idx_end), K.Rank_z(idx_start:idx_end)]; % vector from right ankle -> right ASI hip marker
L_legLen = nan(size(v_Lleg,1),1);
R_legLen = L_legLen;
for i = 1:size(v_Lleg,1)
    L_legLen(i,1) = sqrt(dot(v_Lleg(i,:), v_Lleg(i,:)));
    R_legLen(i,1) = sqrt(dot(v_Rleg(i,:), v_Rleg(i,:)));
end
L_legL = mean(L_legLen,'omitnan');
R_legL = mean(R_legLen,'omitnan');

% Making plots if plt is set to 'Y...' or 'y...'
if strcmp(plt(1),'Y')==1 || strcmp(plt(1),'y')==1
    figure('Name',[fileName ' Heel Height Check'])
    subplot(2,1,1)
    plot(K.time,K.Lhee_z,'b-')
    hold on;
    plot(K.time(idx_start:idx_end),K.Lhee_z(idx_start:idx_end),'r-')
    xlim([0,50]); % only want first 30s of data
    ylim([0,400]); % only care about vertical position between treadmill (0) and 40cm above
    legend('Full Trial','Est. Region','Location','northeast')
    ylabel('Left Heel Height (mm)')

    subplot(2,1,2)
    plot(K.time,K.Rhee_z,'b-')
    hold on;
    plot(K.time(idx_start:idx_end),K.Rhee_z(idx_start:idx_end),'r-')
    xlim([0,50]); % only want first 30s of data
    ylim([0,400]); % only care about vertical position between treadmill (0) and 40cm above
    ylabel('Right Heel Height (mm)')
    xlabel('time (s)')
end