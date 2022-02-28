% This file analyzes the rotational data across speed
plt = 'Y'; % set to 'Y...' or 'y...' to plot check plots
dataDir = 'C:\Users\AJD44\Desktop\Thoracic Rotation\data'; % directory where data files are stored
listing = dir([dataDir, filesep, '*mps_Kin_Filt.csv']); % get listing of files across speeds (only trials with mps in filename)

for i = 1:size(listing,1)
    D_Kin{i,1} = readtable([dataDir filesep listing(i).name]); % Filtered Kinematics from Vicon
    D_Trd{i,1} = readtable([dataDir filesep listing(i).name(1:end-12) 'Trd_Filt.csv']); % Filtered Kinematics from Vicon
    sub_num(i,1) = str2num(listing(i).name(2:3));
    speed{i,1} = movmean(mean([D_Trd{i,1}.Front_vel,D_Trd{i,1}.Rear_vel],2),200)/1000; % smoothed belt speed in meters per second
    [avg_speed{i,1},times{i,1}] = Calc_Avg_Speed(speed{i,1},D_Trd{i,1}.time,listing(i).name,plt);
    [theta_xy{i,1}, theta_ShHip{i,1}, ~, ~] = Calc_Rotation(D_Kin{i,1}, plt, listing(i).name(1:end-13));    
end
