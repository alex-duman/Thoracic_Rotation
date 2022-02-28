function [t_Lstride,t_Rstride,L_strideL,R_strideL] = Calc_Stride_Metrics(D_kin,L_legL,L_legR,plt,trialName,Lheel_height,Rheel_height)
% This function determines the times for each stride cycle as well as the 
% stride length for both the right and left legs, which is defined as when 
% the heel marker reaches it's furthes point forward during the trial.

% Dependencies
% 1) Find_Contact_Events.m (separate function)

% INPUTS
% D_kin - the table of kinematics data associated with that trial
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

% Extracting left and right heel fore-aft motion
L_heel = D_kin.Lhee_y;
R_heel = D_kin.Rhee_y;

T_L = islocalmax(L_heel); % starts of left stride cycles
T_R = islocalmax(R_heel); % starts of right stride cycles

% Inverting data to find local minima
L_heel_neg = -1*D_kin.Lhee_y;
R_heel_neg = -1*D_kin.Rhee_y;
T2_L = islocalmax(L_heel_neg); % furthest distance back of left stride cycles
T2_R = islocalmax(R_heel_neg); % furthest distance back of right stride cycles
t_Lstride_end = D_kin.time(T2_L); % time when left heel is at back of stride cycle
t_Rstride_end = D_kin.time(T2_R); % time when right heel is at back of stride cycle

L_max = L_heel(T_L); % maximum position of left heel in stride cycles
R_max = R_heel(T_R); % maximum position of right heet in stride cycles
L_min = L_heel(T2_L); % minimum position of left heel in stride cycles
R_min = R_heel(T2_R); % minimum position of right heel in stride cycles

% Resize data if needed to only focus on full strides
if length(L_max) ~= length(L_min) % if lengths are not equal
    if length(L_max) > length(L_min) % if more maxima points than minima
        L_max = L_max(1:length(L_min)); % only keep the first maxima associated with minima
    else % more minima points than maxima
        L_min = L_min(1:length(L_max)); % only keep the first minima associated with maxima
    end
end
if length(R_max) ~= length(R_min) % if lengths are not equal
    if length(R_max) > length(R_min) % if more maxima points than minima
        R_max = R_max(1:length(R_min)); % only keep the first maxima associated with minima
    else % more minima points than maxima
        R_min = R_min(1:length(R_max)); % only keep the first minima associated with maxima
    end
end

L_stride_lengths = (L_max-L_min); % distance in millimeters left heel traveled during each stride cycle
R_stride_lengths = (R_max-R_min); % distance in millimeters right heel traveled during each stride cycle
L_strideL = mean(L_stride_lengths)/L_legL; % average left stride length normalized to leg length
R_strideL = mean(R_stride_lengths)/L_legR; % average right stride length normalized to leg length

% Finding the foot touchdown and liftoff times
% LEFT foot
[idx_Lheel_touchdown, idx_Lheel_liftoff] = Find_Contact_Events(D_kin,Lheel_height,plt,'Left',trialName);
% RIGHT foot
[idx_Rheel_touchdown, idx_Rheel_liftoff] = Find_Contact_Events(D_kin,Rheel_height,plt,'Right',trialName);

t_Ltouchdown = D_kin.time(idx_Lheel_touchdown);
t_Lliftoff = D_kin.time(idx_Lheel_liftoff);
t_Rtouchdown = D_kin.time(idx_Rheel_touchdown);
t_Rliftoff = D_kin.time(idx_Rheel_liftoff);

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

if strcmp(plt(1,1),'Y') || strcmp(plt(1,1),'y')
    % PLOTs to check code
    figure('Name',[trialName, ' Stride Metrics'])
    plot(D_kin.time,L_heel,'r-')
    hold on
    plot(D_kin.time,R_heel,'b-')
    plot(D_kin.time(T_L),L_heel(T_L),'rv')
    plot(D_kin.time(T_R),R_heel(T_R),'bv')
    plot(D_kin.time(T2_L),L_heel(T2_L),'r^')
    plot(D_kin.time(T2_R),R_heel(T2_R),'b^')
    legend('left heel','right heel','L heel front','R heel front','L heel back','R heel back','Location','east')
    text(2,mean([R_max; R_min]),{['Left stride = ', num2str(round(L_strideL*L_legL)), ' mm']; ['Right stride = ', num2str(round(R_strideL*L_legR)), ' mm']})
end

end