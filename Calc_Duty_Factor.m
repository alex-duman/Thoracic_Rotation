function [df] = Calc_Duty_Factor(t_Lstride,t_Rstride)
% This function calculates the duty factor for a period based on the
% touchdown and liftoff times of both legs.

% INPUTS
% t_Lstride - table containing the touchdown and liftoff times for the left
%             leg in the first and second columns, respectively
% t_Rstride - table containing the touchdown and liftoff times for the
%             right leg in the first and second columns, respectively

% OUPUT
% df - the duty factor, or fraction of time the foot is in contact with the 
%      ground during the duration of the stride (walking is df >= 0.5 and 
%      running is df < 0.5)  

% Left Leg
if t_Lstride{1,1} < t_Lstride{1,2} % first touchdown occurs prior to first liftoff
    L_stride_dur = diff(t_Lstride{:,1}); % stride duration defined by touchdowns
    L_contact_dur = t_Lstride{1:end-1,2} - t_Lstride{1:end-1,1}; % leave out last step since diff will result in 1 fewer row
    df_L = mean(L_contact_dur./L_stride_dur,'omitnan'); % provides the average duty factor for the left leg
else % first event is liftoff followed by touchdown
    L_stride_dur = diff(t_Lstride{:,1}); % stride duration defined by touchdowns but skip first touchdown 
    L_contact_dur = t_Lstride{2:end,2} - t_Lstride{1:end-1,1}; % subtract off subsequent liftoff since it is shifted below 1 row & make sure same size as diff
    df_L = mean(L_contact_dur./L_stride_dur,'omitnan'); % average duty factor for left leg
end
if df_L == 1 || df_L == 0 % occurs when heel or toe marker drops off
    df_L = nan; % set to nan so it will not affect duty factor calculation
end

% Right Leg
if t_Rstride{1,1} < t_Rstride{1,2} % first touchdown occurs prior to first liftoff
    R_stride_dur = diff(t_Rstride{:,1}); % stride duration defined by touchdowns
    R_contact_dur = t_Rstride{1:end-1,2} - t_Rstride{1:end-1,1}; % leave out last step since diff will result in 1 fewer row
    df_R = mean(R_contact_dur./R_stride_dur,'omitnan'); % provides the average duty factor for the left leg
else % first event is liftoff followed by touchdown
    R_stride_dur = diff(t_Rstride{:,1}); % stride duration defined by touchdowns but skip first touchdown 
    R_contact_dur = t_Rstride{2:end,2} - t_Rstride{1:end-1,1}; % subtract off subsequent liftoff since it is shifted below 1 row & make sure same size as diff
    df_R = mean(R_contact_dur./R_stride_dur,'omitnan'); % average duty factor for left leg
end
if df_R == 1 || df_R == 0 % occurs when heel or toe marker drops off
    df_R = nan; % set to nan so it will not affect duty factor calculation
end

df = mean([df_L,df_R],'omitnan');
end