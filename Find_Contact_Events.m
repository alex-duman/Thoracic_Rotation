function [idx_touchdown, idx_liftoff] = Find_Contact_Events(D_kin,heel_height,plt,side,trialName)
% This function finds the touchdown and liftoff indices and if plt is
% set to anything starting with a 'Y' or 'y' will produce a plot showing
% its estimated liftoff and touchdown events on a graph of heel and toe 
% vertical motion.

% INPUTS
% D_kin - the table of kinematics data associated with that trial
% heel_height - the height of the heel in millimeters (taken from the
%               standing trial) which is used to define the events of
%               interest
% plt - string that contains a 'Y...' or 'y...' if you'd like to plot a
%       summary of the estimated events
% side - a string explaining the heel's side on the body (e.g. 'Left' or
%        'Right')
% trialName - string contianing trial name

% OUTPUTS
% idx_touchdown - logical vector with true conditions corresponding to
%                      foot touchdown event times (determined by when heel 
%                      marker passes through vertical height at standing 
%                      and continues moving downward)
% idx_liftoff - logical vector with true conditions corresponding to
%                    foot liftoff event times (determined by when toe 
%                    marker passes reaches vertical minima just priot to  
%                    the upward foot liftoff)

% Determine column names of interest
if strcmp(side,'Left')==1
    column1 = 'Lhee_z';
    column2 = 'Ltoe_z';
elseif strcmp(side,'Right')==1
    column1 = 'Rhee_z';
    column2 = 'Rtoe_z';
end

% Touchdown Events (defined by heel height)
heel_data = D_kin{:,column1} - heel_height; % subtract off heel height during standing frame
idx_touchdown = zeros(size(heel_data));
prev_touchdown = -0.5; % start with previous touchdown occuring more than 0.5s ahead of any step
for i = 3:length(idx_touchdown)-2
    if sign(heel_data(i)) > sign(heel_data(i+1)) && sign(heel_data(i-2)) > sign(heel_data(i+2)) && (prev_touchdown+0.6) < D_kin.time(i) % moving downward through point i and more than 0.6seconds after previous step
        if sign(heel_data(i)) == 0 % index occured exactly at touchdown
            idx_touchdown(i) = 1; % heel touchdown event
        else
            idx_touchdown(i+1) = 1; % use subsequent step since it will be just after touchdown
        end
        prev_touchdown = D_kin.time(i); % time most recent touchdown calculated occured
    else
        % nothing because we'll keep that row set to zero to indicate no event
    end
end
idx_touchdown = logical(idx_touchdown);
% Liftoff Events (defined by toe height)
toe_data = D_kin{:,column2} - heel_height; % subtract off heel height during standing frame
idx_liftoff = zeros(size(toe_data));
toe_data_neg = toe_data.*-1; % flip to use localmax to find the minima
idx_liftoff = islocalmax(toe_data_neg);
% Ensuring there is only one minima detected within a given timepoint
indices = find(idx_liftoff);
t_liftoffs = D_kin.time(indices);
all_minima = toe_data(indices);
for i = 2:length(indices)-1 % for values where idx_liftoff is true (or logical 1)
    if t_liftoffs(i+1) - t_liftoffs(i) < 0.7 || t_liftoffs(i) - t_liftoffs(i-1) < 0.7 % if adjacent minima are within 0.7 second
        % figuring out which is the actual minima prior to takeoff
        local_idx = (t_liftoffs(i) -0.7) < t_liftoffs & t_liftoffs < (t_liftoffs(i) +0.7);
        step_min = min(toe_data(indices(local_idx)));
        idx_liftoff(indices(i)) = false; % set initially to false
        idx_liftoff(indices(all_minima==step_min)) = true; % set regional min to true (1)
    end
end
if t_liftoffs(2) - t_liftoffs(1) < 0.7 % need to check for final liftoff condition
    local_idx = (t_liftoffs(1) +0.7) < t_liftoffs;
    step_min = min(toe_data(indices(local_idx)));
    idx_liftoff(indices(1)) = false; % set initially to false
    idx_liftoff(indices(all_minima==step_min)) = true; % set regional min to true (1)
end
if t_liftoffs(end) - t_liftoffs(end-1) < 0.7 % need to check for final liftoff condition
    local_idx = (t_liftoffs(end) -0.7) < t_liftoffs;
    step_min = min(toe_data(indices(local_idx)));
    idx_liftoff(indices(end)) = false; % set initially to false
    idx_liftoff(indices(all_minima==step_min)) = true; % set regional min to true (1)
end
for i = 1:length(indices)
    if toe_data(indices(i)) > 0 % if local min occurs above resting heel height for some reason
        idx_liftoff(indices(i)) = false; % don't keep that local minima as a true toe off
    end
end

% Plot to check contact events if specified by plt
if strcmp(plt(1,1),'Y')==1 || strcmp(plt(1,1),'y')==1
    figure('Name',[trialName, ' ', side, ' Contact Data'])
    plot(D_kin.time,D_kin{:,column1},'b-')
    hold on;
    plot(D_kin.time,D_kin{:,column2},'r-')
    plot(D_kin.time(idx_touchdown),D_kin{idx_touchdown,column1},'kv')
    plot(D_kin.time(idx_liftoff),D_kin{idx_liftoff,column2},'k^')
    ylabel('Vertical Position (mm)')
    xlabel('time (s)')
    legend('Heel height','Toe height','touchdown','liftoff','Location','east')
end