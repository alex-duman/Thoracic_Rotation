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
%            each muscle (maximum was taken after the moving average was
%            applied to the signal).

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
max_EMGs = table(max_y(1,1),max_y(1,2),max_y(1,3),'VariableNames',{['RF_',side(1)], ['MG_',side(1)], ['LG_',side(1)]});
percent = [0:0.1:100]';
avg_rf_sup = mean(EMG_rf_sup,2); % getting the average EMG activity for all support phases
avg_mg_sup = mean(EMG_mg_sup,2);
avg_lg_sup = mean(EMG_lg_sup,2);
avg_rf_swi = mean(EMG_rf_swi,2); % getting the average EMG activity for all swing phases
avg_mg_swi = mean(EMG_mg_swi,2);
avg_lg_swi = mean(EMG_lg_swi,2);
varNames = {'percent',['avg_rf_sup_',side(1)], ['avg_mg_sup_',side(1)], ['avg_lg_sup_',side(1)],['avg_rf_swi_',side(1)],['avg_mg_swi_',side(1)],['avg_lg_swi_',side(1)]}; % so body side is included in curve (column) title/name
avg_curves = table(percent,avg_rf_sup,avg_mg_sup,avg_lg_sup,avg_rf_swi,avg_mg_swi,avg_lg_swi,'VariableNames',varNames);

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
        p2(i) = patchline(percent, EMG_rf_swi(:,i),'edgecolo',[0.6,0.6,0.6],'edgealpha',0.2);
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