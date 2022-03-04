function [idx_touchdown, idx_liftoff] = Find_Contact_Events(K,heel_height,plt,side,trialName)
% This function finds the touchdown and liftoff indices and if plt is
% set to anything starting with a 'Y' or 'y' will produce a plot showing
% its estimated liftoff and touchdown events on a graph of heel and toe 
% vertical motion.

% INPUTS
% K - the table of kinematics data associated with that trial
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
heel_data = K{:,column1} - heel_height; % subtract off heel height during standing frame
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
    if sign(heel_data(i)) > sign(heel_data(i+1)) && sign(heel_data(i-2)) > sign(heel_data(i+2)) && heel_data(i-n_frames) > (heel_height/2) % moving downward through point i and n_frames previously heel was at least twice the standing height.
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
% Liftoff Events (defined by toe height)
toe_data = K{:,column2} - heel_height; % subtract off heel height during standing frame
for i = 0:size(touchdowns,1)
    if i == 0 % want minimum prior to first touchdown (may or may not be actual liftoff event)
        [toe_minima(i+1,1), idx_toe_min(i+1,1)] = min(toe_data(1:touchdowns(i+1)));
    elseif i == size(touchdowns,1) % want minimum after last touchdown (may or may not be actual liftoff event)
        [toe_minima(i+1,1), idx_toe_min(i+1,1)] = min(toe_data(touchdowns(i):end));
    else % to find toe minimums between each touchdown (should be actual liftoff events)
        [toe_minima(i+1,1), idx_toe_min(i+1,1)] = min(toe_data(touchdowns(i):touchdowns(i+1)));
    end
end
indices = idx_toe_min+[0;(touchdowns-1)]; % indices for liftoff events, touchdowns are -1 since idx_toe_min = 1 corresponds to exact value of touchdown so we want it to be that same touchdown value
toe_data = [toe_data; zeros(40,1)]; % adding additional 0.4s to end of data at height zero to catch any final touchdowns
for j = [1,size(indices,1)] % want to check the first and last possible liftoff events to ensure they are real
    if toe_minima(j) < 0 && toe_data(indices(j)+40) > 0.1*heel_height % when minimum is below the standing heel height and 0.4s later the foot is 10% above standing condition
        % nothing because this is likely a true liftoff event
    else
        indices(j) = nan; % remove the value since it is not likely to be a true liftoff event
    end
end
indices = rmmissing(indices); % removes any minima that were not actual liftoff events

idx_liftoff = zeros(size(toe_data));
for i = indices
    idx_liftoff(i) = 1; % set index equal to 1 where liftoff events occur
end
idx_liftoff = logical(idx_liftoff);

% Plot to check contact events if specified by plt
if strcmp(plt(1,1),'Y')==1 || strcmp(plt(1,1),'y')==1
    figure('Name',[trialName, ' ', side, ' Contact Data'])
    plot(K.time,K{:,column1},'b-')
    hold on;
    plot(K.time,K{:,column2},'r-')
    plot(K.time(idx_touchdown),K{idx_touchdown,column1},'kv')
    plot(K.time(idx_liftoff),K{idx_liftoff,column2},'k^')
    ylabel('Vertical Position (mm)')
    xlabel('time (s)')
    legend('Heel height','Toe height','touchdown','liftoff','Location','east')
end