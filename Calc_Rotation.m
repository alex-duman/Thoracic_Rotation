function [theta_xy, theta_ShHip, avg_ang_disp_xy, avg_ang_disp_ShHip] = Calc_Rotation(K, plt, trialName)
% This function determines the angle of the shoulders in the transverse
% plane (theta_xy) through time, as well as the angle of the shoulders
% relative to the hips (theta_ShHip), and determines the average angular
% diplacement per stride using local maximums and minimums for both methods
% (avg_ang_disp_xy and avg_ang_disp_ShHip).

% INPUTS
% K - the table containing the kinematic data for the trial of interest

% OUTPUTS
% theta_xy - the angle of the shoulders through time in the transvers (xy-)
%            plane (units of degrees)
% theta_ShHip - the angle of the shoulders through time relative to the
%               hips (units of degrees)
% avg_ang_disp_xy - the average angular displacement of the shoulders in 
%                   the transverse plane across all strides in the trial
%                   (units of degrees)
% avg_ang_disp_ShHip - the average angular displacement of the shoulders
%                      relative to the hips across all strides in the trial
%                      (units of degrees)

% Finding markers of interest
Lsho = [K.Lsho_x, K.Lsho_y, K.Lsho_z]; % left shoulder
Rsho = [K.Rsho_x, K.Rsho_y, K.Rsho_z]; % right shoulder
Lasi = [K.Lasi_x, K.Lasi_y, K.Lasi_z]; % left anterior superior iliac spine
Rasi = [K.Rasi_x, K.Rasi_y, K.Rasi_z]; % right anterior superior iliac spine
time = K.time; % vector of time

% Calculating Theta
for i = 1:size(Lsho,1)
    theta_xy(i,1) = atand((Lsho(i,2)-Rsho(i,2))/(Lsho(i,1)-Rsho(i,1))); % only in transverse (xy-) plane
    Hips_xy(i,1) = atand((Lasi(i,2)-Rasi(i,2))/(Lasi(i,1)-Rasi(i,1))); % hips in transverse (xy-) plane
    theta_ShHip(i,1) = theta_xy(i,1) - Hips_xy(i,1); % shoulder angle relative to hip angle in transverse (xy-) plane   
end

% Determining phases based on local maxima and minima
[peaks,locs_peaks] = findpeaks(theta_xy);
[troughs, ~] = findpeaks(-1*theta_xy);
troughs = -1.*troughs;

avg_peak = mean(peaks(peaks > 0)); % exclude any local maximum that occured below zero degrees (not actual maxima)
avg_trough = mean(troughs(troughs < 0)); % exclude anly local minimum that occured above zero degrees (not actual minima)
avg_ang_disp_xy = avg_peak - avg_trough; % total angular displacement for single stride

for i = 1:(size(peaks,1)-1)
    peaks2(i,1) = max(theta_ShHip(locs_peaks(i,1):locs_peaks((i+1),1),1)); % check for peak between two peaks based on cleaner theta_xy signal
    troughs2(i,1) = min(theta_ShHip(locs_peaks(i,1):locs_peaks((i+1),1),1)); % check for minimums between two peaks based on cleaner theta_xy signal
end
avg_peak2 = mean(peaks2);
avg_trough2 = mean(troughs2);
avg_ang_disp_ShHip = avg_peak2 - avg_trough2; % total angular displacement for single stride

if strcmp(plt(1),'Y')==1 || strcmp(plt(1),'y')==1
    % Plot to check angle of shoulder rotation
    figure('Name', [trialName, ' Shoulder Rotation Check'])
    subplot(2,1,1)
    plot(K.time,theta_xy,'Color',[0.6,0.6,0.6])
    hold on
    yline(avg_peak,'r--','label','Average Peak')
    yline(avg_trough,'b--','label','Average Trough')
    text(K.time(round(length(K.time)/4)), mean([avg_peak, avg_trough]), ['Mean Angular Displacement: ' num2str(round(avg_ang_disp_xy,1)), ' deg'])
    title('Rotation within Transverse (xy-) Plane')
    ylabel('Shoulder Rotation (deg)')

    subplot(2,1,2)
    plot(K.time,theta_ShHip,'Color',[0.6,0.6,0.6])
    hold on
    yline(avg_peak2,'r--','label','Average Peak')
    yline(avg_trough2,'b--','label','Average Trough')
    text(K.time(round(length(K.time)/4)), mean([avg_peak2, avg_trough2]), ['Mean Angular Displacement: ' num2str(round(avg_ang_disp_ShHip,1)), ' deg'])
    title('Rotation relative to Hips in Transverse Plane')
    ylabel('Shoulder Rotation (deg)')
    xlabel('time (s)')
end

end