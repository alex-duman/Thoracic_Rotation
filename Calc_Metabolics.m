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
V_CO2 = mean(VCO2,'omitnan')*1000/BM; % units of ml/kg/min
V_O2 = mean(VO2,'omitnan')*1000/BM; % units of ml/kg/min
RER = mean(rer,'omitnan'); % unitless
HR = mean(hr,'omitnan'); % bpm

if strcmp(plt(1),'Y')==1 || strcmp(plt(1),'y')==1
    % Plot to check metabolic data
    figure('Name',[trialName, ' Metabolics Check'])
    subplot(4,1,1)
    plot(t_Exh(Idx_breaths),M.VO2(Idx_breaths)*1000/BM,'b-o')
    hold on;
    yline(V_O2,'b--','label',['Average V_{O2}: ', num2str(round(V_CO2,1)), ' ml O_2 kg^{-1} min^{-1}'])
    ylabel('VO2 (mL O_2/kg/min)')

    subplot(4,1,2)
    plot(t_Exh(Idx_breaths),M.VCO2(Idx_breaths)*1000/BM,'r-o')
    yline(V_CO2,'r--','label',['Average V_{CO2}: ', num2str(round(V_CO2,1)), ' ml CO_2 kg^{-1} min^{-1}'])
    ylabel('VO2 (mL CO_2/kg/min)')

    subplot(4,1,3)
    plot(t_Exh(Idx_breaths),M.RER(Idx_breaths),'m-o')
    yline(RER,'m--','label',['Average RER: ', num2str(round(RER,2))])
    ylabel('RER (<1 is aerobic)')

    subplot(4,1,4)
    plot(t_Exh(Idx_breaths),M.HR(Idx_breaths),'k-o')
    if isnan(HR)~=1
        yline(HR,'k--','label',['Average HR: ' num2str(round(HR,1)), ' bpm'])
    end
    ylabel('Heart Rate (bpm)')
    xlabel('time (s)')
end

end