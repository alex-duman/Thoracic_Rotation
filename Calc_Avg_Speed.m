function [avg_speed,times] = Calc_Avg_Speed(speed,time,filename,plt)
% This function determines where the constant speed occurs (assumes they
% are ~27 seconds in duration and are spaced 30s apart due to acceleration
% phase.

% INPUTS
% speed - the vector of treadmill speed through time the subject experienced
% time - the time vector coresponding to the treadmill events
% filename - the string contianing the file's full name (e.g. 'S01_0to2mps_Kin_Filt.csv')
% plt - determines whether to plot data to visualize if timepoints were
%       calculated correctly. Will only create plots if set to 'Y...' or 'y...'.

% OUTPUTS
% avg_speed - column vector with corresponding average speeds for each of
%             the periods completed with constant speed.
% times - table containing the start and end times corrersponding to the
%         avg_speed in the same row.

speed = speed - mean(speed(1:5000)); % Remove any offset so stationary is zero m/s

belt_speeds = (str2num(filename(5)) + 0.25):0.25:str2num(filename(8));

idx_speed(1,1) = find((speed >= (belt_speeds(1)-0.025)), 1, 'first') + 1.5/time(2); % 1.5 second following the first time speed increases above 1% below expected speed, here 1/time(2) serves as sampling frequency
idx_speed(1,2) = idx_speed(1,1) + 27/time(2); % 27 seconds after the treadmill gets up to expected speed, here 1/time(2) serves as sampling frequency
for j = 2:size(belt_speeds,2)
    idx_speed(j,1:2) = idx_speed(j-1,1:2) + 30/time(2); % speeds changed every 30s, here 1/time(2) serves as sampling frequency
end

% calculate average speed over period & times when corresponding speed
% starts and ends
for j = 1:size(idx_speed,1)
    avg_speed(j,1) = mean(speed(idx_speed(j,1):idx_speed(j,2)));
    start(j,1) = time(idx_speed(j,1));
    finish(j,1) = time(idx_speed(j,2));
end
times = table(start,finish);

% Plot to visualize times when speed is constant
if strcmp(plt(1),'Y')==1 || strcmp(plt(1),'y')==1
    figure('Name',[filename ' Speeds'])
    plot(time,speed)
    hold on;
    plot(time(idx_speed(:,1)), speed(idx_speed(:,1)),'r>')
    plot(time(idx_speed(:,2)), speed(idx_speed(:,2)),'r<')
end