function [] = Analyze_Indv_Trial(trial,Dir)
% This function analyzes an individual thoracic rotation trial specified by
% the trial name and make plots to help visualize the data from a
% particular trial.

% Dependencies


% INPUTS
% trial - string containing base file name (Ex: 

%% Questions
% 1) Does thoracic rotation affect cost of transport?
%       IVs: amount of thoracic rotation (condition and/or average angular
%            displacement of shoulders)
%       DVs: V_O2 and V_CO2

% 2) Does thoracic rotation reduce peak impact force during locomotion?
%       IVs: amount of thoracic rotation
%       DVs: Peak ground reaction force

% 3) Does thoracic rotation decrease intensity of muscle activation in leg
% muscles?
%       IVs: amount of thoracic rotation
%       DVs: muscle intensity of rectus femoris, MG and LG

% 4) Does preferred amount of rotation change after intervention trials?
%       IVs: whether or not interventions (no rotation and forced rotation)
%            were performed prior to trial
%       DVs: amount of thoracic rotation

%% Perform analyses for specific trial

    % Read in all data for trial
    D_Kin = readtable([Dir filesep trial 'Kin_Filt.csv']); % Filtered Kinematics from Vicon
    D_Met = readtable([Dir filesep trial 'Meta.csv']); % Metabolic Data from Q-Track
    D_Trd = readtable([Dir filesep trial 'Trd_Filt.csv']); % Filtered Treadmill data from Vicon
    D_EMG = readtable([Dir filesep trial 'EMG_Raw.csv']); % Raw EMG data from Vicon

    % Extract Data from trial name
    [subject(i,1), movement{i,1}, RotCond{i,1}] = Extract_trialName_Info(trial{i,1});

    % Calculate amount of thoracic rotation
    [theta_xy{i,1}, theta_ShHip{i,1}, avg_ang_disp_xy{i,1}, avg_ang_disp_ShHip{i,1}] = Calc_Rotation(D_Kin{i,1},trial_plt,trial{i,1});

    % Calculate V_O2 and V_CO2
    [V_O2{i,1},V_CO2{i,1},RER{i,1},HR{i,1}] = Calc_Metabolics(D_Met{i,1},BM(subject(i,1)));

    % Calculate peak Ground Reaction Forces
    [peakGRF(i,1)] = Calc_Forces(D_Trd{i,1},BM(subject(i,1)),trial_plt,trial{i,1});

    % Calculate Stride Metrics
    [t_Lstride{i,1},t_Rstride{i,1},L_strideL(i,1),R_strideL(i,1)] = Calc_Stride_Metrics(D_Kin{i,1},L_legL(subject(i,1)),L_legR(subject(i,1)),trial_plt,trial{i,1},L_heel(subject(i,1)),R_heel(subject(i,1)));

    % Calculate muscle intensities
    [I_Lrf(i,1),I_Lrf_support(i,1),I_Lrf_swing(i,1),I_Lmg(i,1),I_Lmg_support(i,1),I_Lmg_swing(i,1),I_Llg(i,1),I_Llg_support(i,1),I_Llg_swing(i,1), n_steps_L(i,1), avg_curves_L{i,1}, max_EMGs_L{i,1}] = Calc_MuscleInt(D_EMG{i,1},t_Lstride{i,1},'Left',trial_plt,trial{i,1});
    [I_Rrf(i,1),I_Rrf_support(i,1),I_Rrf_swing(i,1),I_Rmg(i,1),I_Rmg_support(i,1),I_Rmg_swing(i,1),I_Rlg(i,1),I_Rlg_support(i,1),I_Rlg_swing(i,1), n_steps_R(i,1), avg_curves_R{i,1}, max_EMGs_R{i,1}] = Calc_MuscleInt(D_EMG{i,1},t_Rstride{i,1},'Right',trial_plt,trial{i,1});

    % Append values to table
    main_table = table();


%% FUNCTIONS

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
    text(K.time(length(K.time)/4), mean([avg_peak, avg_trough]), ['Mean Angular Displacement: ' num2str(round(avg_ang_disp_xy,1)), ' deg'])
    title('Rotation within Transverse (xy-) Plane')
    ylabel('Shoulder Rotation (deg)')

    subplot(2,1,2)
    plot(K.time,theta_ShHip,'Color',[0.6,0.6,0.6])
    hold on
    yline(avg_peak2,'r--','label','Average Peak')
    yline(avg_trough2,'b--','label','Average Trough')
    text(K.time(length(K.time)/4), mean([avg_peak2, avg_trough2]), ['Mean Angular Displacement: ' num2str(round(avg_ang_disp_ShHip,1)), ' deg'])
    title('Rotation relative to Hips in Transverse Plane')
    ylabel('Shoulder Rotation (deg)')
    xlabel('time (s)')
end

end

function [V_O2,V_CO2,RER,HR] = Calc_Metabolics(M,BM,plt,trialName)
% This function will calculate the metabolics including the mass-specific 
% volume of oxygen consumed per unit time (V_O2) and the mass-specific
% volume of carbon dioxide consumed per unit time (V_CO2).

% INPUTS
% M - the table of metabolic data for a given trial
% BM - body mass of subject in kg

% OUTPUTS
% V_O2 - the mass-specific volume of oxygen consumed per unit time
% V_CO2 - the mass-specific volume of oxygen consumed per unit time
% RER - the Respiratory Exchange Ratio which is the ratio of CO2 production
%       to O2 consumption and should be less than or equal to 1 for aerobic
%       exercise.

% finding times associated with each and exhale
Idx_Exh = M.Exh == 1; % index when end of exhales occur
t_Exh = M.time(Idx_Exh); % time at end of exhale

% ASSUME that the last 10 seconds of the trial are no longer at specified
% walking/running condition (delay between turning treadmill off and
% stopping data acquisition on Q-Track).
t_final = M.time(end,1) - 10; % don't want last 10 seconds
t_start = M.time(end,1) - 70; % want 60s total of metabolic data

% Getting breaths of interest, 60 seconds of data from 70 seconds prior to
% end of recording up until 10 seconds prior to end of recording (last 
% minute of locomotion data).
Idx_breaths = t_start <= t_Exh & t_Exh <= t_final;
VCO2 = M.VCO2(Idx_breaths); % units of L/min
VO2 = M.VO2(Idx_breaths); % units of L/min
rer = M.RER(Idx_breaths); % unitless
hr = M.HR(Idx_breaths); % bpm

% Converting to Relative VO2 and VCO2 (ml/kg/min) & taking the average
V_CO2 = mean(VCO2)*1000/BM; % units of ml/kg/min
V_O2 = mean(VO2)*1000/BM; % units of ml/kg/min
RER = mean(rer); % unitless
HR = mean(hr); % bpm

if strcmp(plt(1),'Y')==1 || strcmp(plt(1),'y')==1
    % Plot to check metabolic data
    figure('Name',[trialName, ' Metabolics Check'])
    subplot(4,1,1)
    plot(M.time(Idx_breaths),M.VO2(Idx_breaths)*1000/BM,'b-o')
    hold on;
    yline(V_O2,'b--','label','Average V_{O2}')
    ylabel('VO2 (mL O_2/kg/min)')

    subplot(4,1,2)
    plot(M.time(Idx_breaths),M.VCO2(Idx_breaths)*1000/BM,'r-o')
    yline(V_O2,'r--','label','Average V_{CO2}')
    ylabel('VO2 (mL CO_2/kg/min)')

    subplot(4,1,3)
    plot(M.time(Idx_breaths),M.RER(Idx_breaths),'m-o')
    yline(RER,'m--','label','Average RER')
    ylabel('RER (<1 is aerobic)')

    subplot(4,1,4)
    plot(M.time(Idx_breaths),M.HR(Idx_breaths),'k-o')
    yline(HR,'k--','label','Average HR')
    ylabel('Heart Rate (bpm)')
end

end

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

avg_F_max = mean(F_mag(TF));
peakGRF = avg_F_max/(BM*9.81); % average peak impact force in multiples of body weight

if strcmp(plt(1,1),'Y') || strcmp(plt(1,1),'y')
    % PLOT to check code
    figure('Name',[trialName, ' Average Forces'])
    t = 1:length(F_mag);
    plot(t,F_mag)
    hold on;
    plot(t(TF),F_mag(TF),'r*')
    yline(avg_F_max,'k--','label','Average F max')
    ylabel('Force Magnitude (N)'), xlabel('time (ms)')
end

end

function [I_rf,I_rf_sup, I_rf_swi,I_mg,I_mg_sup, I_mg_swi,I_lg,I_lg_sup, I_lg_swi, n_steps, avg_curves, max_EMGs] = Calc_MuscleInt(E,t_stride,side,plt, trialName)
% This function calculates the average intensity of activation during a 
% stride cycle for the muscles of interest (rectus femoris, medial gastroc 
% and lateral gastroc) from the specified side (Left or Right). It provides
% the intensity for the entire stride as well as the support-phase and
% swing-phase components of the stride cycle separately.

% INPUTS
% E - table of the electromyography data
% t_stride - table with times associated with heel touchdown (1st column) 
%            and heel liftoff (2nd column; see Calc_Stride_Metrics() function for more detail)
% side - string indicating which side of the body is being analyzed ('Left'
%        or 'Right')
% trialName - string containing trial name for figure handle if plt is set
%             to 'y...' or 'Y...'

% OUTPUTS
% I_rf - average rectus femoris intensity during entire stride cycle
% I_rf_sup - average rectus femoris intensity during the support phase
% I_rf_swi - average rectus femoris intensity during the swing phase
% I_mg - average medial gastroc intensity during entire stride cycle
% I_mg_sup - average medial gastroc intensity during support phase
% I_mg_swi - average medial gastroc intensity during swing phase
% I_lg - average lateral gastroc intensity during entire stride cycle
% I_lg_sup - average lateral gastroc intensity during support phase
% I_lg_swi - average lateral gastroc intensity during swing phase
% n_steps - the number of full stride cycles recorded for the leg
% avg_curves - table containing the average EMG signal smoothed via moving
%              mean for all muscles during stance and swing phase
% max_EMGs - table containing the maximum magntidue of the EMG activity for
%            each muscle

if isnan(t_stride{end,1})==1 || isnan(t_stride{end,2})==1 % if not complete final stride
    t_stride(end,:) = []; % want to delete off so only measure complete strides
end

% Determine whether first event is a touchdown or a liftoff
if t_stride{1,1} < t_stride{1,2} % first event is a touchdown
    for i = 1:size(t_stride,1)-1 % finding indices for subsequent stride phases
        idx_sup_start(i) = find(E.time==t_stride{i,1}); % start of support phase
        idx_sup_end(i) = find(E.time==t_stride{i,2}); % end of support (liftoff)
        idx_swi_start(i) = idx_sup_end(i) + 1; % start of swing phase (liftoff)
        idx_swi_end(i) = find(E.time==t_stride{i+1,1}) - 1; % end of swing phase (next touchdown)
        idx_start(i) = find(E.time==t_stride{i,1}); % start of support for stride cycle i
        idx_end(i) = find(E.time==t_stride{i+1,1}) - 1; % end of swing phase for stride cycle i
    end
else % first event is a liftoff
    for i = 1:size(t_stride,1)-1 % finding indices for subsequent stride phases
        idx_swi_start(i) = find(E.time==t_stride{i,2}); % start of swing phase (liftoff)
        idx_swi_end(i) = find(E.time==t_stride{i,1}) - 1; % end of swing phase (next touchdown)
        idx_sup_start(i) = find(E.time==t_stride{i,1}); % start of support phase
        idx_sup_end(i) = find(E.time==t_stride{i+1,2}) - 1; % end of support (liftoff)
        idx_start(i) = find(E.time==t_stride{i,2}); % start of swing for stride cycle i
        idx_end(i) = find(E.time==t_stride{i+1,2}) - 1; % end of support phase for stride cycle i
    end
end

% Sum the rectified EMG signal gives the total activity during phase, which
% we want to normalize to the duration of that phase to account for
% different step durations, therefore we divide the sum of rectified
% intensity by the duration. Assumed the actual sampling frequency was
% 1000Hz since every other time point is zeros for all sensors!
dt = E.time(2) - E.time(1);
dt_sup = E.time(idx_sup_end) - E.time(idx_sup_start); % duration of support phase
dt_swi = E.time(idx_swi_end) - E.time(idx_swi_start); % duration of swing phase
dt_step = E.time(idx_end) - E.time(idx_start); % duration of whole stride cycle
for i = 1:size(t_stride,1)-1
    I_rf_sup_step(i,1) = sum(abs(E{idx_sup_start(i):idx_sup_end(i), [side(1) '_RF']}))*dt / dt_sup(i,1); % support phase
    I_mg_sup_step(i,1) = sum(abs(E{idx_sup_start(i):idx_sup_end(i), [side(1) '_MG']}))*dt / dt_sup(i,1);
    I_lg_sup_step(i,1) = sum(abs(E{idx_sup_start(i):idx_sup_end(i), [side(1) '_LG']}))*dt / dt_sup(i,1);
    I_rf_swi_step(i,1) = sum(abs(E{idx_swi_start(i):idx_swi_end(i), [side(1) '_RF']}))*dt / dt_swi(i,1); % swing phase
    I_mg_swi_step(i,1) = sum(abs(E{idx_swi_start(i):idx_swi_end(i), [side(1) '_MG']}))*dt / dt_swi(i,1);
    I_lg_swi_step(i,1) = sum(abs(E{idx_swi_start(i):idx_swi_end(i), [side(1) '_LG']}))*(1/1000) / dt_swi(i,1);
    I_rf_step(i,1) = sum(abs(E{idx_start(i):idx_end(i), [side(1) '_RF']}))*dt / dt_step(i,1);     % whole stride cycle
    I_mg_step(i,1) = sum(abs(E{idx_start(i):idx_end(i), [side(1) '_MG']}))*dt / dt_step(i,1);
    I_lg_step(i,1) = sum(abs(E{idx_start(i):idx_end(i), [side(1) '_LG']}))*dt / dt_step(i,1);
end

n_steps = size(t_stride,1)-1;

% get average intensity for all steps during trial
I_rf_sup = mean(I_rf_sup_step,'omitnan');
I_mg_sup = mean(I_mg_sup_step,'omitnan');
I_lg_sup = mean(I_lg_sup_step,'omitnan');
I_rf_swi = mean(I_rf_swi_step,'omitnan');
I_mg_swi = mean(I_mg_swi_step,'omitnan');
I_lg_swi = mean(I_lg_swi_step,'omitnan');
I_rf = mean(I_rf_step,'omitnan');
I_mg = mean(I_mg_step,'omitnan');
I_lg = mean(I_lg_step,'omitnan');

% Creating EMG series based on percent of phase rather than time since
% steps may occur over different durations
w = 20; % window size for moving average
% get the rectified & smoothed emg signals for each support and swing phase
% for easier visualization when plotting
for i = 1:size(I_rf_sup_step,1)
    rf_sup{1,i} = movmean(abs(E{idx_sup_start(i):idx_sup_end(i), [side(1) '_RF']}), w); % getting smoothed moving average of rectified EMG
    EMG_rf_sup(:,i) = interp1(linspace(0,100,size(rf_sup{1,i},1)), rf_sup{i}, [0:0.1:100]); % interpolating to get 1000 points between 0% and 100% support phase
    mg_sup{i} = movmean(abs(E{idx_sup_start(i):idx_sup_end(i), [side(1) '_MG']}), w);
    EMG_mg_sup(:,i) = interp1(linspace(0,100,size(mg_sup{i},1)), mg_sup{i}, [0:0.1:100]);
    lg_sup{i} = movmean(abs(E{idx_sup_start(i):idx_sup_end(i), [side(1) '_LG']}), w);
    EMG_lg_sup(:,i) = interp1(linspace(0,100,size(lg_sup{i},1)), lg_sup{i}, [0:0.1:100]);
end
for i = 1:size(I_rf_swi_step,1)
    rf_swi{i} = movmean(abs(E{idx_swi_start(i):idx_swi_end(i), [side(1) '_RF']}), w);
    EMG_rf_swi(:,i) = interp1(linspace(0,100,size(rf_swi{i},1)), rf_swi{i}, [0:0.1:100]);
    mg_swi{i} = movmean(abs(E{idx_swi_start(i):idx_swi_end(i), [side(1) '_MG']}), w);
    EMG_mg_swi(:,i) = interp1(linspace(0,100,size(mg_swi{i},1)), mg_swi{i}, [0:0.1:100]);
    lg_swi{i} = movmean(abs(E{idx_swi_start(i):idx_swi_end(i), [side(1) '_LG']}), w);
    EMG_lg_swi(:,i) = interp1(linspace(0,100,size(lg_swi{i},1)), lg_swi{i}, [0:0.1:100]);
end
max_y = [max(max([EMG_rf_sup, EMG_rf_swi])), max(max([EMG_mg_sup, EMG_mg_swi])), max(max([EMG_lg_sup, EMG_lg_swi]))];
max_EMGs = table(max_y(1,1),max_y(1,2),max_y(1,3),'VariableNames',{'RF', 'MG', 'LG'});
percent = 0:0.1:100;
avg_rf_sup = mean(EMG_rf_sup,2); % getting the average EMG activity for all support phases
avg_mg_sup = mean(EMG_mg_sup,2);
avg_lg_sup = mean(EMG_lg_sup,2);
avg_rf_swi = mean(EMG_rf_swi,2); % getting the average EMG activity for all swing phases
avg_mg_swi = mean(EMG_mg_swi,2);
avg_lg_swi = mean(EMG_lg_swi,2);
avg_curves = table(avg_rf_sup,avg_mg_sup,avg_lg_sup,avg_rf_swi,avg_mg_swi,avg_lg_swi);

% PLOT to check average EMG
if strcmp(plt(1),'Y')==1 || strcmp(plt(1),'y')==1
    % Creating Figure
    figure('Name', [trialName, ' ', side, ' EMG Signals'])
    subplot(3,2,1)
    p1(1) = patchline(percent, EMG_rf_sup(:,1),'edgecolor',[0.6,0.6,0.6],'edgealpha',0.2);
    hold on;
    for i =2:size(EMG_rf_sup,2)
        p1(i) = patchline(percent, EMG_rf_sup(:,i),'edgecolo',[0.6,0.6,0.6],'edgealpha',0.2);
    end
    plot(percent,avg_rf_sup,'r-')
    ylim([0 max_y(1,1)]);
    ylabel('RF intensity')
    title('Support Phase')

    subplot(3,2,2)
    p2(1) = patchline(percent, EMG_rf_swi(:,1),'edgecolor',[0.6,0.6,0.6],'edgealpha',0.2);
    hold on;
    for i =2:size(t_stride,1)-1
        p2(i) = patchline(percent, EMG_rf_sup(:,i),'edgecolo',[0.6,0.6,0.6],'edgealpha',0.2);
    end
    ylim([0 max_y(1,1)]);
    plot(percent,avg_rf_swi,'r-')
    title('Swing Phase')
    ylabel('RF')

    subplot(3,2,3)
    p3(1) = patchline(percent, EMG_mg_sup(:,1),'edgecolor',[0.6,0.6,0.6],'edgealpha',0.2);
    hold on;
    for i =2:size(t_stride,1)-1
        p3(i) = patchline(percent, EMG_mg_sup(:,i),'edgecolo',[0.6,0.6,0.6],'edgealpha',0.2);
    end
    ylim([0 max_y(1,2)]);
    plot(percent,avg_mg_sup,'r-')
    ylabel('MG intensity')

    subplot(3,2,4)
    p4(1) = patchline(percent, EMG_mg_swi(:,1),'edgecolor',[0.6,0.6,0.6],'edgealpha',0.2);
    hold on;
    for i =2:size(t_stride,1)-1
        p4(i) = patchline(percent, EMG_mg_swi(:,i),'edgecolo',[0.6,0.6,0.6],'edgealpha',0.2);
    end
    ylim([0 max_y(1,2)]);
    plot(percent,avg_mg_swi,'r-')
    ylabel('MG')

    subplot(3,2,5)
    p5(1) = patchline(percent, EMG_lg_sup(:,1),'edgecolor',[0.6,0.6,0.6],'edgealpha',0.2);
    hold on;
    for i =2:size(t_stride,1)-1
        p5(i) = patchline(percent, EMG_lg_sup(:,i),'edgecolo',[0.6,0.6,0.6],'edgealpha',0.2);
    end
    ylim([0 max_y(1,3)]);
    plot(percent,avg_lg_sup,'r-')
    ylabel('LG intensity')
    xlabel('Support Phase (%)')

    subplot(3,2,6)
    p6(1) = patchline(percent, EMG_lg_swi(:,1),'edgecolor',[0.6,0.6,0.6],'edgealpha',0.2);
    hold on;
    for i =2:size(t_stride,1)-1
        p6(i) = patchline(percent, EMG_lg_swi(:,i),'edgecolo',[0.6,0.6,0.6],'edgealpha',0.2);
    end
    ylim([0 max_y(1,3)]);
    plot(percent,avg_lg_swi,'r-')
    ylabel('LG')
    xlabel('Swing Phase (%)')
end

end

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

function [subject, movement, RotCond] = Extract_trialName_Info(trialName)
% This function extracts the information we want to use in order to compare
% the data, like the subject number, movement condition and rotational
% condition.

% INPUTS
% trialName - the string containing the trial name of interest (Ex: S01_wlk_Preferred_pre_)

% OUTPUTS
% subject - the assigned subject number for the participant of that trial (Ex: 1)
% movement - the movement classification (Ex: 'walking')
% RotCond - rotational condition (Ex: 'Preferred Pre')

% Subject number
subject = str2num(trialName(1,2:3));

% Movement Condition
if strcmp(trialName(1,5:7),'wlk') % if trial is labeled as walking
    movement = 'walk';
else % not labeled as walking
    movement = nan;
    warning(['Trial ' trialName(1,1:end-1) ' does not have a recognized condition!']);
end

% Rotational Condition
if strcmp(trialName(1,9:10),'No') % No rotation condition
    RotCond = 'None';
elseif strcmp(trialName(1,9:14),'Forced') % Forced rotation condition
    RotCond = 'Forced';
elseif strcmp(trialName(1,end-4:end-1),'pre') % Preferred Pre
    RotCond = 'Preferred Pre';
elseif strcmp(trialName(1,end-4:end-1),'pos') % Preferred Post
    RotCond = 'Preferred Post';
else % none of the above conditions
    RotCond = nan;
    warning(['Trial ' trialName(1,1:end-1) ' does not have a recognized rotational condition!']);
end

end