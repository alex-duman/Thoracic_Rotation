function [peakGRF] = Calc_Forces(F,BM,plt,trialName)
% This function will calculate the average peak ground reaction force
% (peakGRF) for a given trial.

% INPUTS
% F - the table of treadmill data containing the force plate measurements
% BM - body mass of subject in kg
% plt - whether or not to plot the force data (only plots if set to 'Y',
% 'Yes', 'y', or 'yes')
% trialName - name of the trial to use for figure title

% OUTPUTS
% peakGRF - the average peak impact force during locomotion across all
%             steps of the trial in multiples of Body Weight

Force = [F.Front_Fx, F.Front_Fy, F.Front_Fz];

% getting magnitude of impact force by using matrix algebra to take square 
% root of the sum of the individual components at each time point (row of F)
F_mag = sqrt(diag(Force*Force')); % magnitude of impact force from front (impacting) force plate

TF = islocalmax(F_mag) & F_mag > BM*9.81; % Find local maximums and only select those that occur when force is greater than 1 body weight
indices = find(TF);

% Removing any false positive local maxima (should not occur within 0.25s)
for i = 1:size(indices,1)-1
    if F.time(indices(i+1)) - F.time(indices(i)) < 0.5 % want to only keep local maxima that are true max forces from a step not from double contact phases
        if F_mag(indices(i)) < F_mag(indices(i+1))
            TF(indices(i)) = false;
        else % indices(i) is the true local maxima for the step
            TF(indices(i+1)) = false;
        end
    end
end

avg_F_max = mean(F_mag(TF));
peakGRF = avg_F_max/(BM*9.81); % average peak impact force in multiples of body weight

if strcmp(plt(1,1),'Y') || strcmp(plt(1,1),'y')
    % PLOT to check code
    figure('Name',[trialName, ' Average Forces'])
    t = 1:length(F_mag);
    plot(t,F_mag)
    hold on;
    plot(t(TF),F_mag(TF),'r*')
    yline(avg_F_max,'k--','label',['Average F max: ' num2str(round(avg_F_max)), ' N.'])
    ylabel('Force Magnitude (N)'), xlabel('time (ms)')
end

end