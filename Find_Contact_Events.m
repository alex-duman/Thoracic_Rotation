function [idx_touchdown, idx_liftoff] = Find_Contact_Events(K,heel_height,toe_height,plt,side,trialName)
% This function finds the touchdown and liftoff indices and if plt is
% set to anything starting with a 'Y' or 'y' will produce a plot showing
% its estimated liftoff and touchdown events on a graph of heel and toe 
% vertical motion.

% INPUTS
% K - the table of kinematics data associated with that trial
% heel_height - the height of the heel in millimeters (taken from the
%               standing trial) which is used to define the events of
%               interest
% toe_height - the height of the toe in millimeters (taken from standing
%              trial) which is used to define takeoff events, as well as
%              landing events for toe-runners.
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
heel_data = K{:,column1} - 1.1*heel_height; % subtract off 110% heel height during standing frame (touchdown occurs when heel height is ~110% of standing height)
% add zeros at beginning of heel_data to help capture any
% steps occuring in first few timepoints measured (assumed filming
% frequency of 100 Hz)
if str2num(trialName(5:8)) < 0.4 % for slow trials check up to 50 frames before since they are slow moving
    n_frames = 50;
else % quick enough where 30 frames (0.3s) is sufficient
    n_frames = 30;
end
heel_data = [zeros(n_frames,1); heel_data];
idx_touchdown = zeros(size(heel_data));
for i = (n_frames+3):length(idx_touchdown)-2
    if sign(heel_data(i)) > sign(heel_data(i+1)) && sign(heel_data(i-2)) > sign(heel_data(i+2)) && heel_data(i-n_frames) > (heel_height/2) % moving downward through point i and n_frames previously heel was at 1.5x's the standing height.
        if sum(idx_touchdown) >= 1 && (i-find(idx_touchdown,1,'last')) > n_frames % when touchdown event already found, make sure subsequent touchdown is at least n_frames after previous
            if sign(heel_data(i)) == 0 % index occured exactly at touchdown
                idx_touchdown(i) = 1; % heel touchdown event
            else
                idx_touchdown(i+1) = 1; % use subsequent step since it will be just after touchdown
            end
        elseif sum(idx_touchdown) == 0 % first touchdown
            if sign(heel_data(i)) == 0 % index occured exactly at touchdown
                idx_touchdown(i) = 1; % heel touchdown event
            else
                idx_touchdown(i+1) = 1; % use subsequent step since it will be just after touchdown
            end
        end
    else
        % nothing because we'll keep that row set to zero to indicate no event
    end
end
idx_touchdown = idx_touchdown((n_frames+1):end); % remove first n_frames that were artificially added
idx_touchdown = logical(idx_touchdown);
touchdowns = find(idx_touchdown);
TD_id = 'heel';

% IF fewer than 4 touchdown events then it's likely TOE-Runner data rather
% than a heel-striker
if size(touchdowns,1) < 4 % || max(diff(touchdowns)) > 1.5*mean(diff(touchdowns))
% recognized less than 4 steps using heel to identify touchdowns (poor metric and likely toe-runner)
% OR if the maximum difference between touchdown is greater then 1.5x's the
% average difference (likely missing one or more steps and toe may offer
% better results.
    toe_dataTD = K{:,column2} - toe_height; % subtract off toe height during standing frame
    toe_dataTD = [zeros(n_frames,1); toe_dataTD];
    idx_touchdown = zeros(size(toe_dataTD));
    for i = (n_frames+3):length(idx_touchdown)-2
        if sign(toe_dataTD(i)) > sign(toe_dataTD(i+1)) && sign(toe_dataTD(i-2)) > sign(toe_dataTD(i+2)) && toe_dataTD(i-n_frames) > 10 % moving downward through point i and n_frames previously toe was at least 10mm above the standing height.
            if sum(idx_touchdown) >= 1 && (i-find(idx_touchdown,1,'last')) > n_frames % when touchdown event already found, make sure subsequent touchdown is at least n_frames after previous
                if sign(toe_dataTD(i)) == 0 % index occured exactly at touchdown
                    idx_touchdown(i) = 1; % toe touchdown event
                else
                    idx_touchdown(i+1) = 1; % use subsequent step since it will be just after touchdown
                end
            elseif sum(idx_touchdown) == 0 % first touchdown
                if sign(toe_dataTD(i)) == 0 % index occured exactly at touchdown
                    idx_touchdown(i) = 1; % toe touchdown event
                else
                    idx_touchdown(i+1) = 1; % use subsequent step since it will be just after touchdown
                end
            end
        else
            % nothing because we'll keep that row set to zero to indicate no event
        end
    end
    idx_touchdown = idx_touchdown((n_frames+1):end); % remove first n_frames that were artificially added
    idx_touchdown = logical(idx_touchdown);
    touchdowns = find(idx_touchdown);
    TD_id = 'toe';
end

% Liftoff Events (defined by toe height)
toe_data = K{:,column2} - toe_height; % subtract off toe height during standing frame (threshold for contact is when toe is 100% above standing contact)
liftoffs = nan(size(touchdowns,1)+1,1);
LO_id = 'toe';
for i = 0:size(touchdowns,1)
    if i == 0
        toe = toe_data(1:touchdowns(1),1); % get toe data prior to first touchdown (trial may start with takeoff)
    elseif i == size(touchdowns,1)
        toe = toe_data(touchdowns(end):end,1); % get toe data after last touchdown (trial may end with takeoff)
    else
        toe = toe_data(touchdowns(i):touchdowns(i+1),1); % get toe data between each touchdown
    end

    % Determining when takeoff event occurs
    for j = 11:(size(toe,1)-10)
        if sign(toe(j,1)) < sign(toe(j+1,1)) && sign(toe(j-10,1)) < sign(toe(j+10,1)) && toe(j+10,1) > 0 % moving upward through point j and 10 frames after toe is still above standing height.
            if i == 0
                liftoffs(i+1,1) = j;
            else
                liftoffs(i+1,1) = j + touchdowns(i);
                TD_LO_difference(i,1) = j; % frames between touchdown and subsequent liftoff
            end
        end
    end
end
if sum(isnan(liftoffs)) >= size(liftoffs,1)/2 % if more than half of liftoff events are missing (likely toe marker dropped off)
    LO_id = 'heel';
    heel_dataLO = K{:,column1} - heel_height; % subtract off toe height during standing frame
    for i = 0:size(touchdowns,1)
        if i == 0
            heel = heel_dataLO(1:touchdowns(1),1); % get heel data prior to first touchdown (trial may start with takeoff)
        elseif i == size(touchdowns,1)
            heel = heel_dataLO(touchdowns(end):end,1); % get heel data after last touchdown (trial may end with takeoff)
        else
            heel = heel_dataLO(touchdowns(i):touchdowns(i+1),1); % get toe data between each touchdown
        end

        % Determining when takeoff event occurs
        for j = 6:(size(heel,1)-5)
            if sign(heel(j,1)) < sign(heel(j+1,1)) && sign(heel(j-5,1)) < sign(heel(j+5,1)) && heel(j+5,1) > 0 && heel(j,1) < heel_height % moving upward through point j and 5 frames after toe is still above standing height.
                if i == 0
                    liftoffs(i+1,1) = j;
                else
                    liftoffs(i+1,1) = j + touchdowns(i);
                    TD_LO_difference(i,1) = j; % frames between touchdown and subsequent liftoff
                end
            end
        end
    end
end
TD_LO_avg_diff = round(mean(TD_LO_difference,'omitnan')); % average amount of frames between touchdown and subsequent toeoff
for i = 2:size(liftoffs,1)-1
    if isnan(liftoffs(i,1))==1 % If could not find value
        liftoffs(i,1) = touchdowns(i-1,1) + TD_LO_avg_diff; % Estimates takeoff based on average of all other steps if could not identify
    end
end
if isnan(liftoffs(1,1))==1 % No liftoff event prior to first touchdown
    liftoffs = liftoffs(2:end,1); % remove first potential liftoff since it is NaN
end
if isnan(liftoffs(end,1))==1 % No liftoff event after last touchdown
    liftoffs = liftoffs(1:end-1,1); % remove last potential liftoff since it is NaN
end
idx_liftoff = zeros(size(toe_data));
for i = 1:size(liftoffs)
    idx_liftoff(liftoffs(i),1) = 1;
end
idx_liftoff = logical(idx_liftoff);

% Plot to check contact events if specified by plt
if strcmp(plt(1,1),'Y')==1 || strcmp(plt(1,1),'y')==1
    figure('Name',[trialName, ' ', side, ' Contact Data'])
    plot(K.time,K{:,column1},'b-')
    hold on;
    plot(K.time,K{:,column2},'r-')
    if strcmp(TD_id,'heel')==1 % heel was used for touchdown identification
        plot(K.time(idx_touchdown),K{idx_touchdown,column1},'kv')
    else % toe was used for touchdown identification
        plot(K.time(idx_touchdown),K{idx_touchdown,column2},'kv')
    end
    if strcmp(LO_id,'toe')==1 % toe used for liftoff identification
        plot(K.time(idx_liftoff),K{idx_liftoff,column2},'k^')
    else % heel used for liftoff identification
        plot(K.time(idx_liftoff),K{idx_liftoff,column1},'k^')
    end
    ylabel('Vertical Position (mm)')
    xlabel('time (s)')
    legend('Heel height','Toe height','touchdown','liftoff','Location','east')
end