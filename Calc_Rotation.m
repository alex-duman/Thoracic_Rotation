function [theta_xy, theta_ShHip, avg_ang_disp_xy, avg_ang_disp_ShHip] = Calc_Rotation(K, tLstride, tRstride, plt, trialName)
% This function determines the angle of the shoulders in the transverse
% plane (theta_xy) through time, as well as the angle of the shoulders
% relative to the hips (theta_ShHip), and determines the average angular
% diplacement per stride using local maximums and minimums for both methods
% (avg_ang_disp_xy and avg_ang_disp_ShHip).

% INPUTS
% K - the table containing the kinematic data for the trial of interest
% tLstride - table containing the touchdown times (1st column) and liftoff
%             times (2nd column) for each stride cycle (rows) of the left leg
% tRstride - table containing the touchdown times (1st column) and liftoff
%             times (2nd column) for each stride cycle (rows) of the left leg
% plt - string to determine whether or not to visualize data through a
%       plot, set to 'Y...' or 'y...' if you want to plot
% trialName - string of full trial name to use for figure name

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

% Calculating Theta (orientation of shoulders)
for i = 1:size(Lsho,1)
    theta_xy(i,1) = atand((Lsho(i,2)-Rsho(i,2))/(Lsho(i,1)-Rsho(i,1))); % only in transverse (xy-) plane
    Hips_xy(i,1) = atand((Lasi(i,2)-Rasi(i,2))/(Lasi(i,1)-Rasi(i,1))); % hips in transverse (xy-) plane
    theta_ShHip(i,1) = theta_xy(i,1) - Hips_xy(i,1); % shoulder angle relative to hip angle in transverse (xy-) plane   
end

%% Determining phases based on strides

% Want to use the foot that starts a full cycle first (greater potential
% for more strides)
if tLstride{1,1} < tRstride{1,1}
    t_stride = tLstride;
    option = 1; % maximum rotational values will occur near touchdown
else % right leg starts a stride cycle prior to left
    t_stride = tRstride;
    option = 2; % minimum rotational values will occur near touchdown
end

if option == 1 % Left leg used and touchdown occur near maximal rotational values
    % Maximums
    for i = 1:size(t_stride,1)
        if ~isnan(t_stride{i,1})
            idx_near_max(i,1) = find(K.time==t_stride{i,1});
            if i == 1
                [theta_xy_max(i,1), idx_xy_max(i,1)] = max(theta_xy(1:(idx_near_max(i,1)+20)));
                [theta_ShHip_max(i,1), idx_ShHip_max(i,1)] = max(theta_ShHip(1:(idx_near_max(i,1))+20)); % and check 0.2 seconds following
            elseif i == size(t_stride,1)
                [theta_xy_max(i,1), idx_xy_max(i,1)] = max(theta_xy((idx_near_max(i,1)-20):end));
                [theta_ShHip_max(i,1), idx_ShHip_max(i,1)] = max(theta_ShHip((idx_near_max(i,1)-20):end));% and check 0.2 second before
                idx_xy_max(i,1) = idx_xy_max(i,1) - 20; % remove 20 from scanning over additional 0.2 seconds
                idx_ShHip_max(i,1) = idx_ShHip_max(i,1) - 20; % remove 20 from scanning over additional 0.2 seconds
            else
                half_idx = round((idx_near_max(i,1)-idx_near_max(i-1,1))/2); % half the index distance between strides
                [theta_xy_max(i,1), idx_xy_max(i,1)] = max(theta_xy(idx_near_max(i-1,1)+half_idx:idx_near_max(i,1)+half_idx));
                [theta_ShHip_max(i,1), idx_ShHip_max(i,1)] = max(theta_ShHip(idx_near_max(i-1,1)+half_idx:idx_near_max(i,1)+half_idx));
                idx_xy_max(i,1) = idx_xy_max(i,1) - half_idx; % subtract off addtional indices added prior to idx_near_max
                idx_ShHip_max(i,1) = idx_ShHip_max(i,1) - half_idx; % subtract off addtional indices added prior to idx_near_max
            end
        else
            idx_near_max(i,1) = nan;
            theta_xy_max(i,1) = nan;
            idx_xy_max(i,1) = nan;
            theta_ShHip_max(i,1) = nan;
            idx_ShHip_max(i,1) = nan;
        end
    end
    % Minimums
    for i = 1:(size(t_stride,1)-1)
        idx_start(i,1) = find(K.time==t_stride{i,1});
        if ~isnan(t_stride{i+1,1})
            idx_end(i,1) = find(K.time==t_stride{i+1,1});
            [theta_xy_min(i,1), idx_xy_min(i,1)] = min(theta_xy(idx_start(i,1):idx_end(i,1)));
            [theta_ShHip_min(i,1), idx_ShHip_min(i,1)] = min(theta_ShHip(idx_start(i,1):idx_end(i,1)));
        else
            idx_end(i,1) = nan;
            theta_xy_min(i,1) = nan;
            idx_xy_min(i,1) = nan;
            theta_ShHip_min(i,1) = nan;
            idx_ShHip_min(i,1) = nan;
        end
    end
    idx_Rotxy_max = [0; idx_near_max(2:end)] + idx_xy_max;
    idx_Rotxy_min = idx_start + idx_xy_min;
    idx_RotSH_max = [0; idx_near_max(2:end)] + idx_ShHip_max;
    idx_RotSH_min = idx_start + idx_ShHip_min;
else % option 2: Right leg used and touchdown occur near minimal rotational values
    % Minimums
    for i = 1:size(t_stride,1)
        if ~isnan(t_stride{i,1})
            idx_near_min(i,1) = find(K.time==t_stride{i,1});
            if i == 1
                [theta_xy_min(i,1), idx_xy_min(i,1)] = min(theta_xy(1:(idx_near_min(i,1)+20)));
                [theta_ShHip_min(i,1), idx_ShHip_min(i,1)] = min(theta_ShHip(1:(idx_near_min(i,1))+20)); % and check 0.2 seconds following
            elseif i == size(t_stride,1)
                [theta_xy_min(i,1), idx_xy_min(i,1)] = min(theta_xy((idx_near_min(i,1)-20):end));
                [theta_ShHip_min(i,1), idx_ShHip_min(i,1)] = min(theta_ShHip((idx_near_min(i,1)-20):end));% and check 0.2 second before
                idx_xy_min(i,1) = idx_xy_min(i,1) - 20; % remove 20 from scanning over additional 0.2 seconds
                idx_ShHip_min(i,1) = idx_ShHip_min(i,1) - 20; % remove 20 from scanning over additional 0.2 seconds
            else
                half_idx = round((idx_near_min(i,1)-idx_near_min(i-1,1))/2); % half the index distance between strides
                [theta_xy_min(i,1), idx_xy_min(i,1)] = min(theta_xy(idx_near_min(i-1,1)+half_idx:idx_near_min(i,1)+half_idx));
                [theta_ShHip_min(i,1), idx_ShHip_min(i,1)] = min(theta_ShHip(idx_near_min(i-1,1)+half_idx:idx_near_min(i,1)+half_idx));
                idx_xy_min(i,1) = idx_xy_min(i,1) - half_idx; % subtract off addtional indices added prior to idx_near_max
                idx_ShHip_min(i,1) = idx_ShHip_min(i,1) - half_idx; % subtract off addtional indices added prior to idx_near_max
            end
        else
            idx_near_min(i,1) = nan;
            theta_xy_min(i,1) = nan;
            idx_xy_min(i,1) = nan;
            theta_ShHip_min(i,1) = nan;
            idx_ShHip_min(i,1) = nan;
        end
    end
    % Maximums
    for i = 1:(size(t_stride,1)-1)
        idx_start(i,1) = find(K.time==t_stride{i,1});
        if ~isnan(t_stride{i+1,1})
            idx_end(i,1) = find(K.time==t_stride{i+1,1});
            [theta_xy_max(i,1), idx_xy_max(i,1)] = max(theta_xy(idx_start(i,1):idx_end(i,1)));
            [theta_ShHip_max(i,1), idx_ShHip_max(i,1)] = max(theta_ShHip(idx_start(i,1):idx_end(i,1)));
        else
            idx_end(i,1) = nan;
            theta_xy_max(i,1) = nan;
            idx_xy_max(i,1) = nan;
            theta_xy_ShHip_max(i,1) = nan;
            idx_ShHip_max(i,1) = nan;
        end
    end
    idx_Rotxy_max = idx_start + idx_xy_max;
    idx_Rotxy_min = [0; idx_near_min(2:end)] + idx_xy_min;
    idx_RotSH_max = idx_start + idx_ShHip_max;
    idx_RotSH_min = [0; idx_near_min(2:end)] + idx_ShHip_min;
end
max_length_xy = min([size(theta_xy_max,1),size(theta_xy_min,1)]); % find the max length for addition/subtraction of vectors
avg_ang_disp_xy = mean(theta_xy_max(1:max_length_xy)-theta_xy_min(1:max_length_xy),'omitnan'); % this way only subsequent rotations are subtracted
max_length_SH = min([size(theta_ShHip_max,1),size(theta_ShHip_min,1)]); % find the max length for addition/subtraction of vectors
avg_ang_disp_ShHip = mean(theta_ShHip_max(1:max_length_SH)-theta_ShHip_min(1:max_length_SH),'omitnan');

%% Visualizing with Plots (if plt set to 'Y...' or 'y...')
if strcmp(plt(1),'Y')==1 || strcmp(plt(1),'y')==1
    % Plot to check angle of shoulder rotation
    figure('Name', [trialName, ' Shoulder Rotation Check'])
    subplot(2,1,1)
    plot(K.time,theta_xy,'Color',[0.6,0.6,0.6])
    hold on
    for i = 1:max_length_xy
        if ~isnan(idx_Rotxy_max(i,1)) && ~isnan(idx_Rotxy_min(i,1))
            plot(K.time(idx_Rotxy_max(i,1)), theta_xy(idx_Rotxy_max(i,1)),'rv')
            plot(K.time(idx_Rotxy_min(i,1)), theta_xy(idx_Rotxy_min(i,1)),'r^')
        end
    end
    text(K.time(round(length(K.time)/4)), 0, ['Mean Angular Displacement: ' num2str(round(avg_ang_disp_xy,1)), ' deg'])
    title('Rotation within Transverse (xy-) Plane')
    ylabel('Shoulder Rotation (deg)')

    subplot(2,1,2)
    plot(K.time,theta_ShHip,'Color',[0.6,0.6,0.6])
    hold on
    for i = 1:max_length_SH
        if ~isnan(idx_RotSH_max(i,1)) && ~isnan(idx_RotSH_min(i,1))
            plot(K.time(idx_RotSH_max(i,1)), theta_ShHip(idx_RotSH_max(i,1)),'rv')
            plot(K.time(idx_RotSH_min(i,1)), theta_ShHip(idx_RotSH_min(i,1)),'r^')
        end
    end
    text(K.time(round(length(K.time)/4)), 0, ['Mean Angular Displacement: ' num2str(round(avg_ang_disp_ShHip,1)), ' deg'])
    title('Rotation relative to Hips in Transverse Plane')
    ylabel('Shoulder Rotation (deg)')
    xlabel('time (s)')
end

end