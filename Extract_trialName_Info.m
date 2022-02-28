function [subject, movement, RotCond] = Extract_trialName_Info(trialName)
% This function extracts the information we want to use in order to compare
% the data, like the subject number, movement condition and rotational
% condition.

% INPUTS
% trialName - the string containing the trial name of interest (Ex: S01_wlk_Preferred_pre_)

% OUTPUTS
% subject - the assigned subject number for the participant of that trial (Ex: 1)
% movement - the movement classification string inside cell (Ex: {'walking'})
% RotCond - rotational condition string inside cell (Ex: {'Preferred Pre'})

% Subject number
subject = str2num(trialName(1,2:3));

% Movement Condition
if strcmp(trialName(1,5:7),'wlk') % if trial is labeled as walking
    movement = {'walk'};
else % not labeled as walking
    movement = nan;
    warning(['Trial ' trialName(1,1:end-1) ' does not have a recognized condition!']);
end

% Rotational Condition
if strcmp(trialName(1,9:10),'No') % No rotation condition
    RotCond = {'None'};
elseif strcmp(trialName(1,9:14),'Forced') % Forced rotation condition
    RotCond = {'Forced'};
elseif strcmp(trialName(1,end-3:end-1),'pre') % Preferred Pre
    RotCond = {'Preferred Pre'};
elseif strcmp(trialName(1,end-3:end-1),'pos') % Preferred Post
    RotCond = {'Preferred Post'};
else % none of the above conditions
    RotCond = nan;
    warning(['Trial ' trialName(1,1:end-1) ' does not have a recognized rotational condition!']);
end

end