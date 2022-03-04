% This script allows you to analyze all the data, from every trial from 
% each individual that has files in the data folder.

trial_plt = 'y'; % set to 'Y...' or 'y...' to see data visualization checks for individual trials
Indv_plt = 'y'; % set to 'Y...' or 'y...' to see data visualization for all trials from one subject together
final_plt = 'y'; % set to 'Y...' or 'y...' to get overall plots comparing all individuals across the different conditions
statistics = 'y'; % set to 'Y...' or 'y...' when you want to analyze the statistics of the outputed RESULTS table
opts = detectImportOptions('C:\Users\AJD44\Desktop\Thoracic Rotation\data\Subject_Morphometrics.csv','NumHeaderLines',2);
opts.VariableNamesLine = 1; opts.VariableUnitsLine = 2; opts.DataLine = 3;
M = readtable('C:\Users\AJD44\Desktop\Thoracic Rotation\data\Subject_Morphometrics.csv',opts); % Morphometrics Table
dataDir = 'C:\Users\AJD44\Desktop\Thoracic Rotation\data'; % directory where data files are stored

for i = 1:size(M,1)
    subject = M.Subject(i);
    if exist('RESULTS','var') == 0
        RESULTS = Analyze_All_Trials_for_Indv(subject,dataDir,M,trial_plt,Indv_plt);
    else
        New_result = Analyze_All_Trials_for_Indv(subject,dataDir,M,trial_plt,Indv_plt);
        RESULTS = [RESULTS; New_result];
    end
end

disp('Finished Analyzing All Data!');

if strcmp(final_plt(1),'Y') || strcmp(final_plt(1),'y')
% Making Plots
RESULTS.RotCond = categorical(RESULTS.RotCond);
RESULTS.RotCond = reordercats(RESULTS.RotCond,{'None'; 'Forced'; 'Preferred Pre'; 'Preferred Post'});

ylabels = {'Rotation (deg)';'Rotation (deg)';'VO_2 (ml O_2/kg/min)';'VCO_2 (ml CO_2/kg/min)';'RER';'Heart Rate (bpm)';'Peak Impact Force (BW)';'Left Stride Length (leg length)';'Right Stride Length (leg length)'};

for i = 4:12
    figure('Name',['Summary ', RESULTS.Properties.VariableNames{i}])
    boxplot(RESULTS{:,i},RESULTS.RotCond,'positions',[0.75:1:4],'Widths',0.3,'Symbol','*r')
    hold on;
    hAx = gca; % retrieve the axes handle
    xtk = hAx.XTick; % retrieve XTicks
    for j = unique(RESULTS.subject)'
        Indv_Results = sortrows(RESULTS(RESULTS.subject==j,:),'RotCond'); % sorted to be in order of xTicks (order of categories)
        plot(xtk +0.25,Indv_Results{:,i},'-','Color',[0.8,0.8,0.8])
        plot(xtk +0.25,Indv_Results{:,i},'k.','MarkerSize',16)
    end
    title(RESULTS.Properties.VariableNames{i})
    ylabel(ylabels{i-3,1});
    xlabel('Thoracic Rotation Condition')
end
end

if strcmp(statistics(1),'Y') || strcmp(statistics(1),'y')
    modelStrings = { ' ~ 1 + (1|subject)' ... % model containing only the subject as random effect
                 ' ~ 1 + RotCond + (1|subject)'}; % model with only 1st order effect

statsVars = RESULTS.Properties.VariableNames([4:12]); % only compares variables of interest in columns 4-12

logLkhdMatrix = nan(length(modelStrings)-1,length(statsVars));

for j = 1:length(statsVars) % loop through all measured variables to test
    cVname = statsVars{j};
    for k = 1:length(modelStrings) %loop through  models to compare for hypothesis testing
        cModelStr = modelStrings{k};
        compLME{k,j} = fitlme(RESULTS,[cVname cModelStr]);    %#ok<*AGROW>
        AICmatrix(k,j) = compLME{k,j}.ModelCriterion{1, 1};
        if k>1
            t_loglklyhood = compare(compLME{k-1,j}, compLME{k,j});
            logLkhdMatrix(k-1,j) = t_loglklyhood.pValue(2);   
        end
    end
end

[bestModelAIC,bestModelIDx] = min(AICmatrix,[],1); %find model with lowest AIC

% Get stats for the best model
modelToUse = mode(bestModelIDx); %find model that modal lowest AIC across variables 

for j = 1:length(statsVars)
    stats{1,j} = anova(compLME{modelToUse,j});
    rSquared(j) = compLME{modelToUse,j}.Rsquared;
    compLME{modelToUse,j} %#ok<NOPTS>
    disp(statsVars{j});
    disp(stats{1,j})
    [b_rand{1,j},bnames{1,j},randStats{1,j}] = randomEffects(compLME{modelToUse,j});
end

for j = 1:length(statsVars)
    cVname = statsVars{j};   
    % Calculate LME model for the current variables   
    t_coeff = compLME{modelToUse,j}.Coefficients;       
    coefficients{1,j} = t_coeff;   
    % Get the random effect results
    % Remove the effect of individual to run posthoc comparison of the fixed effects alone
    yFit_c = fitted(compLME{modelToUse,j},'Conditional',false); % Fit including fixed and random effects
    yFit_m = fitted(compLME{modelToUse,j},'Conditional',true); % Contribution from only fixed effects (marginal)
    y_ind = yFit_m - yFit_c; % Take the difference between the two above to get the effect of individual only. 
    
    % Calculate the corrected y-variable,  with the effect of individual removed
    cData_res = RESULTS.(cVname) - y_ind; 
    
    % Run an ANOVA model with 'anovan' to get the outputs necessary 
    % for the posthoc pairwise comparisons
   [pVals_ANOVAN{j},anovaTable{j},anovaStats{j}, anovaTerms{j}] = anovan(cData_res,...
       {RESULTS.RotCond}...
        ,'random',[], 'continuous',[], 'model', 'interaction', 'varnames',...
        {'condNum'}, 'display', 'off');
    [pairwiseComps_c,MeansSE_c,~,groupNames_c] = multcompare(anovaStats{j},'ctype', 'bonferroni','Dimension',[1],'display','off');
    pairWisePvals{j} = [pairwiseComps_c(:,6)];
end
end