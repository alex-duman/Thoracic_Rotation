function [Indv_table] = Analyze_All_Trials_for_Indv(subject,dataDir,M,trial_plt,Indv_plt)
% This function will run all the trials for a single individual and will
% make summary plots for the individual comparing the different conditions
% when prompted via Indv_plt.

% INPUTS
% subject - integer associated with subject number in study
% dataDir - directory path where data is stored
% M - table of morphometric data for each subject
% trial_plt - if set to 'Y...' or 'y...' it will make summary plots to
%             visualize the data within the trial
% Indv_plt - if set to 'Y...' or 'y...' it will make summary plots to
%             visualize the average data for each condition for that subject

% OUTPUTS
% Indv_table - table containing all data of interest for the individual of
%              interst

% Getting list of files associated with particular subject
SubNum = num2str(subject);
if subject < 10 % need to add zero before single digit subject numbers
    SubNum = ['0', SubNum];
end
if strcmp(M.Movement{1,1}(1),'w')==1
    gait = '_wlk'; % level walking
elseif strcmp(M.Movement{1,1}(1),'i')==1
    gait = '_inc'; % inclined walking
elseif strcmp(M.Movement{1,1}(1),'r')==1
    gait = '_run'; % running
else
    warning('Gait not recognized for ', num2str(M.Subject(1)), ' ', M.Movement{1,1},'!');
end

listing = dir([dataDir, filesep,'S', SubNum, gait, '*Meta.csv']); % get listing of files that have metabolic data (only have one for each trial)

row = 1; % setup initial row of final subject output table
for i = 1:size(listing,1) % number of total trials for all subjects
    % Read in all data for trial
    trial{i,1} = listing(i).name(1,1:end-8);
    [Indv_table(i,:),EMG_curves{1,i},max_EMGs(i,:),RotAngle_curves{1,i}] = Analyze_Trial(trial{i,1},dataDir,M,trial_plt);
end

if strcmp(Indv_plt(1),'Y')==1 || strcmp(Indv_plt(1),'y')==1
    line_style = {'r-';'b-';'k-';'k--'};
    % Thoracic Rotation of Shoulders
    figure('Name',['Subject', SubNum, ' Thoracic Rotation Comparison'])
    subplot(2,1,1)
    plot(RotAngle_curves{1,4}.time,RotAngle_curves{1,4}.theta_xy,line_style{4}); % since alphabetical last one should be Preferred_pre condition
    hold on;
    for i = 1:size(listing,1)-1 % don't need to include last trial since this is Preferred_pre we plotted initially
        plot(RotAngle_curves{1,i}.time,RotAngle_curves{1,i}.theta_xy,line_style{i});
    end
    ylabel('Shoulder Rotation (deg)')
    title('Rotation within Transverse (xy-) Plane')
    legend('Preferred Pre','Forced Rotation','No Rotation','Preferred Post')

    subplot(2,1,2)
    plot(RotAngle_curves{1,4}.time,RotAngle_curves{1,4}.theta_ShHip,line_style{4}); % since alphabetical last one should be Preferred_pre condition
    hold on;
    for i = 1:size(listing,1)-1 % don't need to include last trial since this is Preferred_pre we plotted initially
        plot(RotAngle_curves{1,i}.time,RotAngle_curves{1,i}.theta_ShHip,line_style{i});
    end
    ylabel('Shoulder Rotation (deg)')
    xlabel('time (s)')
    title('Rotation relative to Hips in Transverse Plane')

    % EMG
    max_EMG = max(max_EMGs{:,4:end});
    for i = 1:size(EMG_curves,2)
        for m = 0:2 % number of muscles (rf, mg, lg)
            EMG_curves{1,i}{:,([2,5]+m)} = EMG_curves{1,i}{:,([2,5]+m)}./max_EMG(1+m); % normalize to max intensity recorded on Left
            EMG_curves{1,i}{:,([8,11]+m)} = EMG_curves{1,i}{:,([8,11]+m)}./max_EMG(4+m); % normalize to max intensity recorded on Right
        end
    end
    % Want to find max value for normalized EMG recording so concatenating
    % tables and finding max value for each muscle in each leg
    for i=1:size(EMG_curves,2)
        if i == 1
            All_EMG_curves = EMG_curves{1,i};
        else
            All_EMG_curves = [All_EMG_curves; EMG_curves{1,i}];
        end
    end
    for m = 0:2 % number of muscles (rf, mg, lg)
        max_y(1,m+1) = max(max(All_EMG_curves{:,[2,5]+m})); % left
        max_y(2,m+1) = max(max(All_EMG_curves{:,[8,11]+m})); % right
    end
    Legs = {'Left','Right'};
    for j = 1:size(Legs,2)
    figure('Name',['Subject', SubNum, ' ', Legs{1,j},' EMG Comparison'])
    subplot(3,2,1)
    plot(EMG_curves{1,4}.percent, EMG_curves{1,4}{:,['avg_rf_sup_', Legs{1,j}(1)]},line_style{4})
    hold on;
    for i = 1:size(listing,1)-1 % don't need to include last trial since this is Preferred_pre we plotted initially
        plot(EMG_curves{1,i}.percent, EMG_curves{1,i}{:,['avg_rf_sup_', Legs{1,j}(1)]},line_style{i})
    end
    ylabel('RF activity')
    ylim([0,max_y(j,1)])
    title('Support Phase')

    subplot(3,2,2)
    plot(EMG_curves{1,4}.percent, EMG_curves{1,4}{:,['avg_rf_swi_', Legs{1,j}(1)]},line_style{4})
    hold on;
    for i = 1:size(listing,1)-1 % don't need to include last trial since this is Preferred_pre we plotted initially
        plot(EMG_curves{1,i}.percent, EMG_curves{1,i}{:,['avg_rf_swi_', Legs{1,j}(1)]},line_style{i})
    end
    ylabel('RF activity')
    ylim([0,max_y(j,1)])
    title('Swing Phase')

    subplot(3,2,3)
    plot(EMG_curves{1,4}.percent, EMG_curves{1,4}{:,['avg_mg_sup_', Legs{1,j}(1)]},line_style{4})
    hold on;
    for i = 1:size(listing,1)-1 % don't need to include last trial since this is Preferred_pre we plotted initially
        plot(EMG_curves{1,i}.percent, EMG_curves{1,i}{:,['avg_mg_sup_', Legs{1,j}(1)]},line_style{i})
    end
    ylabel('MG activity')
    ylim([0,max_y(j,2)])
    
    subplot(3,2,4)
    plot(EMG_curves{1,4}.percent, EMG_curves{1,4}{:,['avg_mg_swi_', Legs{1,j}(1)]},line_style{4})
    hold on;
    for i = 1:size(listing,1)-1 % don't need to include last trial since this is Preferred_pre we plotted initially
        plot(EMG_curves{1,i}.percent, EMG_curves{1,i}{:,['avg_mg_swi_', Legs{1,j}(1)]},line_style{i})
    end
    ylabel('MG activity')
    ylim([0,max_y(j,2)])
    legend('Preferred Pre','Forced Rotation','No Rotation','Preferred Post','Location','best')
    
        subplot(3,2,5)
    plot(EMG_curves{1,4}.percent, EMG_curves{1,4}{:,['avg_lg_sup_', Legs{1,j}(1)]},line_style{4})
    hold on;
    for i = 1:size(listing,1)-1 % don't need to include last trial since this is Preferred_pre we plotted initially
        plot(EMG_curves{1,i}.percent, EMG_curves{1,i}{:,['avg_lg_sup_', Legs{1,j}(1)]},line_style{i})
    end
    ylabel('LG activity')
    ylim([0,max_y(j,3)])
    xlabel('Support Phase (%)')

    subplot(3,2,6)
    plot(EMG_curves{1,4}.percent, EMG_curves{1,4}{:,['avg_lg_swi_', Legs{1,j}(1)]},line_style{4})
    hold on;
    for i = 1:size(listing,1)-1 % don't need to include last trial since this is Preferred_pre we plotted initially
        plot(EMG_curves{1,i}.percent, EMG_curves{1,i}{:,['avg_lg_swi_', Legs{1,j}(1)]},line_style{i})
    end
    ylabel('LG activity')
    ylim([0,max_y(j,3)])
    xlabel('Swing Phase (%)')
end

end