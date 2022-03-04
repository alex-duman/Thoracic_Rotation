function [trial_table,EMG_curves,max_EMGs,RotAngle_curves] = Analyze_Trial(trial,dataDir,M,trial_plt)
% This function analyzes an individual thoracic rotation trial specified by
% the trial name and make plots to help visualize the data from a
% particular trial.

% Dependencies


% INPUTS
% trial - string containing base file name (Ex: 'S01_wlk_ForcedRot_001_')
% dataDir - directory where data files are stored
% M - table containing morphometric data from all subjects where each row
%     is the same number associated with the subjects' numbering
% trial_plt - if set to 'Y...' or 'y...' it will make summary plots to
%             visualize the data within the trial

% OUTPUTS
% trial_table - table containing the main variables of interest
% EMG_curves - table containing the average EMG curves for all steps in the 
%              muscles of both legs with a 'percent' column relating to the
%              percent of support or swing phase.
% max_EMGs - table containing the maximum magntidue of the EMG activity for
%            each muscle observed over all steps in the trial (maximum was 
%            taken after the moving average was applied to the signal).

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
    D_Kin = readtable([dataDir filesep trial 'Kin_Filt.csv']); % Filtered Kinematics from Vicon
    D_Met = readtable([dataDir filesep trial 'Meta.csv']); % Metabolic Data from Q-Track
    D_Trd = readtable([dataDir filesep trial 'Trd_Filt.csv']); % Filtered Treadmill data from Vicon
    D_EMG = readtable([dataDir filesep trial 'EMG_Raw.csv']); % Raw EMG data from Vicon

    % Extract Data from trial name
    [subject, movement, RotCond] = Extract_trialName_Info(trial);

    % Calculate Stride Metrics
    [t_Lstride,t_Rstride,L_strideLen,R_strideLen,~] = Calc_Stride_Metrics(D_Kin,M.L_legL(subject),M.L_legR(subject),trial_plt,trial,M.L_heel(subject),M.R_heel(subject));

    % Calculate amount of thoracic rotation
    [theta_xy, theta_ShHip, avg_ang_disp_xy, avg_ang_disp_ShHip] = Calc_Rotation(D_Kin, t_Lstride, t_Rstride, trial_plt, trial);

    % Calculate V_O2 and V_CO2
    [V_O2,V_CO2,RER,HR] = Calc_Metabolics(D_Met,M.BM(subject),trial_plt,trial);

    % Calculate peak Ground Reaction Forces
    [peakGRF] = Calc_Forces(D_Trd,M.BM(subject),trial_plt,trial);

    % Calculate muscle intensities
    [I_Lrf,I_Lrf_support,I_Lrf_swing,I_Lmg,I_Lmg_support,I_Lmg_swing,I_Llg,I_Llg_support,I_Llg_swing, n_steps_L, avg_curves_L, max_EMGs_L] = Calc_MuscleInt(D_EMG,t_Lstride,'Left',trial_plt,trial);
    [I_Rrf,I_Rrf_support,I_Rrf_swing,I_Rmg,I_Rmg_support,I_Rmg_swing,I_Rlg,I_Rlg_support,I_Rlg_swing, n_steps_R, avg_curves_R, max_EMGs_R] = Calc_MuscleInt(D_EMG,t_Rstride,'Right',trial_plt,trial);

    % Append values to table
    trial_table = table(subject,movement,RotCond,avg_ang_disp_xy,avg_ang_disp_ShHip,V_O2,V_CO2,RER,HR,peakGRF,L_strideLen,R_strideLen,I_Lrf,I_Lrf_support,I_Lrf_swing,I_Lmg,I_Lmg_support,I_Lmg_swing,I_Llg,I_Llg_support,I_Llg_swing, n_steps_L,I_Rrf,I_Rrf_support,I_Rrf_swing,I_Rmg,I_Rmg_support,I_Rmg_swing,I_Rlg,I_Rlg_support,I_Rlg_swing, n_steps_R);
    % Average Time-series traces to output for comparison against other trials
    RotAngle_curves = table(D_Kin.time, theta_xy, theta_ShHip,'VariableNames',{'time', 'theta_xy','theta_ShHip'});
    EMG_curves = [avg_curves_L, removevars(avg_curves_R,{'percent'})];
    max_EMGs = [table(subject,movement,RotCond), max_EMGs_L, max_EMGs_R];

end