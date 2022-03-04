function [t_Lstride,t_Rstride,L_strideL,R_strideL,df] = Calc_Stride_Metrics(K,L_legL,L_legR,plt,trialName,Lheel_height,Rheel_height)
% This function determines the times for each stride cycle as well as the 
% stride length for both the right and left legs, which is defined as when 
% the heel marker reaches it's furthes point forward during the trial.

% Dependencies
% 1) Find_Contact_Events.m (separate function)

% INPUTS
% K - the table of kinematics data associated with that trial
% L_legL - left leg length in millimeters
% L_legR - right leg length in millimeters
% plt - whether or not to plot the force data (only plots if set to 'Y',
% 'Yes', 'y', or 'yes')
% trialName - name of the trial to use for figure title
% Lheel_height - the subject's left heel height in millimeters (taken from
%                standing trial)
% Rheel_height - the subject's right heel height in millimeters

% OUTPUTS
% L_strideL - average left stride length
% R_strideL - average right stride length
% t_Lstride - table containing the touchdown times (1st column) and liftoff
%             times (2nd column) for each stride cycle (rows) of the left leg
% t_Rstride - table containing the touchdown times (1st column) and liftoff
%             times (2nd column) for each stride cycle (rows) of the left leg
% df - duty factor (walking is df >= 0.5 and running is df < 0.5)

% Finding the foot touchdown and liftoff times
% LEFT foot
[idx_Lheel_touchdown, idx_Lheel_liftoff] = Find_Contact_Events(K,Lheel_height,plt,'Left',trialName);
% RIGHT foot
[idx_Rheel_touchdown, idx_Rheel_liftoff] = Find_Contact_Events(K,Rheel_height,plt,'Right',trialName);

t_Ltouchdown = K.time(idx_Lheel_touchdown);
t_Lliftoff = K.time(idx_Lheel_liftoff);
t_Rtouchdown = K.time(idx_Rheel_touchdown);
t_Rliftoff = K.time(idx_Rheel_liftoff);

if length(t_Ltouchdown) == length(t_Lliftoff)
    t_Lstride = table(t_Ltouchdown,t_Lliftoff);
elseif length(t_Ltouchdown) > length(t_Lliftoff)
    t_Lliftoff = [t_Lliftoff; nan]; % add nan to end of liftoffs to get to same length
    t_Lstride = table(t_Ltouchdown,t_Lliftoff);
elseif length(t_Ltouchdown) < length(t_Lliftoff)
    t_Ltouchdown = [t_Ltouchdown; nan]; % add nan to end of liftoffs to get to same length
    t_Lstride = table(t_Ltouchdown,t_Lliftoff);
end
if length(t_Rtouchdown) == length(t_Rliftoff)
    t_Rstride = table(t_Rtouchdown,t_Rliftoff);
elseif length(t_Rtouchdown) > length(t_Rliftoff)
    t_Rliftoff = [t_Rliftoff; nan]; % add nan to end of liftoffs to get to same length
    t_Rstride = table(t_Rtouchdown,t_Rliftoff);
elseif length(t_Rtouchdown) < length(t_Rliftoff)
    t_Rtouchdown = [t_Rtouchdown; nan]; % add nan to end of liftoffs to get to same length
    t_Rstride = table(t_Rtouchdown,t_Rliftoff);
end

% Calculate Duty Factor
df = Calc_Duty_Factor(t_Lstride,t_Rstride);

% Calculating Stide Lengths
[L_strideL, L_idx_max, L_idx_min] = Calc_Stride_Length(K, 'Left', t_Lstride, L_legL);
[R_strideL, R_idx_max, R_idx_min] = Calc_Stride_Length(K, 'Right', t_Rstride, L_legR);

if strcmp(plt(1,1),'Y') || strcmp(plt(1,1),'y')
    % PLOTs to check code
    figure('Name',[trialName, ' Stride Metrics'])
    plot(K.time,K.Lhee_y,'r-')
    hold on
    plot(K.time,K.Rhee_y,'b-')
    plot(K.time(L_idx_max(1,1)),K.Lhee_y(L_idx_max(1,1)),'rv') % plot i = 1 first for legend to be in order
    plot(K.time(L_idx_min(1,1)),K.Lhee_y(L_idx_min(1,1)),'r^')
    plot(K.time(R_idx_max(1,1)),K.Rhee_y(R_idx_max(1,1)),'bv')
    plot(K.time(R_idx_min(1,1)),K.Rhee_y(R_idx_min(1,1)),'b^')
    for i = 2:size(L_idx_max,1)
        if ~isnan(L_idx_max(i,1))
            plot(K.time(L_idx_max(i,1)),K.Lhee_y(L_idx_max(i,1)),'rv')
        end
    end
    for i = 2:size(L_idx_min,1)
        if ~isnan(L_idx_min(i,1))
            plot(K.time(L_idx_min(i,1)),K.Lhee_y(L_idx_min(i,1)),'r^')
        end
    end
    for i = 2:size(R_idx_max,1)
        if ~isnan(R_idx_max(i,1))
            plot(K.time(R_idx_max(i,1)),K.Rhee_y(R_idx_max(i,1)),'bv')
        end
    end
    for i = 2:size(R_idx_min,1)
        if ~isnan(R_idx_min(i,1))
            plot(K.time(R_idx_min(i,1)),K.Rhee_y(R_idx_min(i,1)),'b^')
        end
    end
    legend('left heel','right heel','L heel front','L heel back','R heel front','R heel back','Location','east')
    text(K.time(1),mean(K.Lhee_y,'omitnan'),{['Left stride = ', num2str(round(L_strideL*L_legL)), ' mm']; ['Right stride = ', num2str(round(R_strideL*L_legR)), ' mm']})
end

end