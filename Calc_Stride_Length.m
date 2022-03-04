function [strideLength, idx_max, idx_min] = Calc_Stride_Length(K, side, t_stride, legLen)
% This function determines the stride length for a given leg based on the
% maximum and minimal positions of the heel in the fore-aft direciton
% (y-axis for Vicon).

% INPUTS
% K - table of kinematics data from Vicon
% side - string containing either 'R...' or 'L...' for the leg side being
%        analyzed
% t_stride - table of times corresponding to touchdown and liftoff times
%            for that leg in the first and second columns, respectively
% legLen - length of leg in millimeters

% OUTPUTS
% strideLength - average stride length for all strides measured
% idx_max - indeces for stride maximal fore-aft positions
% idx_min - indices for stride minimal fore-aft positions

% Extracting heel fore-aft motion
heel = K{:,[side(1), 'hee_y']};
QuartStep = round(mean(abs(t_stride{:,1} - t_stride{:,2}),'omitnan')/4*100); % average # of frames for 1/4 step so can find max and min within 1/2 step area

% Stride Lengths
for i = 1:size(t_stride,1)
    if ~isnan(t_stride{i,1}) % touchdown events will be close to maximal fore-aft position of heel during stride
        idx_max(i,1) = find(K.time==t_stride{i,1}); % touchdown event i
        if idx_max(i,1)-QuartStep > 0 && idx_max(i,1)+QuartStep <= size(heel,1) % have data quarter step in either direction
            [Step_max(i,1),idx_max_mod] = max(heel(idx_max(i,1)-QuartStep:idx_max(i,1)+QuartStep)); % maximum fore-aft heel position occurs near touchdown (generally slightly after)
            idx_max(i,1) = idx_max(i,1) + (idx_max_mod-(QuartStep+1)); % provides actual index of maximum
        elseif idx_max(i,1)-QuartStep < 0 % don't have enough data before
            [Step_max(i,1),idx_max(i,1)] = max(heel(1:idx_max(i,1)+QuartStep)); % check for max from beginning
        else % don't have enough data afterward
            [Step_max(i,1),idx_max_mod] = max(heel(idx_max(i,1)-QuartStep:end)); % check for max up until end of speed condition
            idx_max(i,1) = idx_max(i,1) + (idx_max_mod-(QuartStep+1)); % provides actual index of maximum
        end
    else % not a complete stride
        idx_max(i,1) = nan;
        Step_max(i,1) = nan;
    end
    if ~isnan(t_stride{i,2}) % liftoff events will be the minimum fore-aft position of foot during stride
        idx_min(i,1) = find(K.time==t_stride{i,2}); % liftoff event i
        if idx_min(i,1)-QuartStep > 0 && idx_min(i,1)+QuartStep <= size(heel,1) % have data quarter step in either direction
            [Step_min(i,1), idx_min_mod] = min(heel(idx_min(i,1)-QuartStep:idx_min(i,1)+QuartStep)); % minimum fore-aft heel position occurs near liftoff (usually slightly after)
            idx_min(i,1) = idx_min(i,1) + (idx_min_mod-(QuartStep+1)); % provides actual index of minimum
        elseif idx_min(i,1)-QuartStep < 0 % don't have enough data before
            [Step_min(i,1), idx_min(i,1)] = min(heel(1:idx_min(i,1)+QuartStep)); % check for minimum from beginning
        else % don't have enough data afterward
            [Step_min(i,1), idx_min_mod] = min(heel(idx_min(i,1)-QuartStep:end)); % check for minimum until end of speed condition
            idx_min(i,1) = idx_min(i,1) + (idx_min_mod-(QuartStep+1)); % provides actual index of minimum
        end
    else
        idx_max(i,1) = nan;
        Step_min(i,1) = nan;
    end
end
strideLength = mean([Step_max - Step_min],'omitnan')/legLen; % average left stride length normalized to leg length

end