% This file analyzes the rotational data across speed
plt = 'N'; % set to 'Y...' or 'y...' to plot check plots
dataDir = 'C:\Users\AJD44\Desktop\Thoracic Rotation\data'; % directory where data files are stored
listing = dir([dataDir, filesep, '*mps_Kin_Filt.csv']); % get listing of files across speeds (only trials with mps in filename)

for i = 1:size(listing,1)
    D_Kin{i,1} = readtable([dataDir filesep listing(i).name]); % Filtered Kinematics from Vicon
    D_Trd{i,1} = readtable([dataDir filesep listing(i).name(1:end-12) 'Trd_Filt.csv']); % Filtered Kinematics from Vicon
    sub_num(i,1) = str2num(listing(i).name(2:3));
    speed{i,1} = movmean(mean([D_Trd{i,1}.Front_vel,D_Trd{i,1}.Rear_vel],2),200)/1000; % smoothed belt speed in meters per second
    [avg_speed{i,1},times{i,1}] = Calc_Avg_Speed(speed{i,1},D_Trd{i,1}.time,listing(i).name(1:end-13),plt);
    % Get stride data to better determine thoracic rotation
    [Lheel_height(i,1),Rheel_height(i,1),Ltoe_height(i,1), Rtoe_height(i,1),L_legL(i,1),R_legL(i,1)] = Calc_Foot_Height(D_Kin{i,1},speed{i,1},D_Trd{i,1}.time,plt,listing(i).name(1:end-13));
    for j = 1:size(times{i,1},1)
        trialName = [listing(i).name(1:4), num2str(round(avg_speed{i,1}(j),2)), ' mps'];
        idx{i,j} = (D_Kin{i,1}.time >= times{i,1}{j,1}) & (D_Kin{i,1}.time <= times{i,1}{j,2});
        [t_Lstride{i,j},t_Rstride{i,j},L_strideL{i,j},R_strideL{i,j},df{i,1}(j,1)] = Calc_Stride_Metrics(D_Kin{i,1}(idx{i,j},:),L_legL(i,1),R_legL(i,1),plt,trialName,Lheel_height(i,1),Rheel_height(i,1),Ltoe_height(i,1),Rtoe_height(i,1));
        [theta_xy{i,j}, theta_ShHip{i,j}, avg_Rot_xy{i,1}(j,1), avg_Rot_ShHip{i,1}(j,1)] = Calc_Rotation(D_Kin{i,1}(idx{i,j},:), t_Lstride{i,j}, t_Rstride{i,j}, plt, trialName);
    end
end

subjects = unique(sub_num);
for j = 1:size(subjects,1)
    for i = 1:size(listing,1)
        sub_idx = sub_num == subjects(j);
        dutyfactor{j,1} = cat(1, df{sub_idx,1});
        sub_speeds{j,1} = cat(1, avg_speed{sub_idx,1});
        Rotation{j,1} = [cat(1,avg_Rot_xy{sub_idx,1}), cat(1,avg_Rot_ShHip{sub_idx,1})]; % Rotation of shoulders in xy-plane in 1st column and relative to hips in xy-plane in 2nd column
    end
end

colors = [127,0,0; 179,0,0; 215,48,31; 239,101,72; 252,141,89; 
          253,187,132; 253,212,158; 254,232,200; 255,247,236; 255,255,255;
          0,0,127; 0,0,179; 31,48,215; 72,101,239; 89,141,252;
          132,187,253; 158,212,253; 200,232,254; 236,247,255; 0,0,0]./255;


% Plotting
figure('Name','Thoracic Rotation Across Speed')
subplot(2,2,1)
plot(sub_speeds{1,1}(1,1),Rotation{1,1}(1,1),'o','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',colors(1,:))
hold on;
for i = 1:size(subjects,1)
    for j = 1:size(Rotation{i,1},1)
        if dutyfactor{i,1}(j,1) >= 0.5 % walking
            plot(sub_speeds{i,1}(j,1),Rotation{i,1}(j,1),'o','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',colors(i,:))
        else % running
            plot(sub_speeds{i,1}(j,1),Rotation{i,1}(j,1),'s','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',colors(i,:))
        end
    end
    for k = 1:2 % splitting walking vs running
        if k == 1 % walking
            index = dutyfactor{i,1} >= 0.5;
        else % running
            index = dutyfactor{i,1} < 0.5;
        end
        % Linear Regresssion (assuming intercept should be Rot = 0)
        B_xy(:,i,k) = [ones(size(sub_speeds{i,1}(index,1))), sub_speeds{i,1}(index,1)]\Rotation{i,1}(index,1); % from Y = X[B0; B1], here B0 is the y-intercept and B1 is the slope
        y_reg_xy{i,k} = B_xy(2,i,k)*sub_speeds{i,1}(index,1) + B_xy(1,i,k); % y-values for the linear regression line
        plot(sub_speeds{i,1}(index,1),y_reg_xy{i,k},'-','Color',colors(i,:))
    end
end
ylabel({'Shoulder Rotation';'xy-plane (deg)'})

% Determine average and standard deviation for walking and running in
% xy-plane
X = categorical({'Walking','Running'});
X = reordercats(X,{'Walking','Running'});
wlk_slopes_xy = B_xy(2,:,1);
run_slopes_xy = B_xy(2,:,2);
Y_slope_xy_avg = [mean(wlk_slopes_xy,'omitnan'), mean(run_slopes_xy,'omitnan')];
subplot(2,2,2)
bar(X,Y_slope_xy_avg)
hold on;
er = errorbar(X,Y_slope_xy_avg,[std(wlk_slopes_xy,'omitnan'),std(run_slopes_xy,'omitnan')]);
er.Color = [0,0,0]; er.LineStyle = 'none';
X_rand_wlk = 0.1*(rand(size(wlk_slopes_xy))-0.5) + 1;
X_rand_run = 0.1*(rand(size(run_slopes_xy))-0.5) + 2;
plot(X_rand_wlk,wlk_slopes_xy,'.','Color',[1, 0.3, 0.3])
plot(X_rand_run,run_slopes_xy,'.','Color',[1, 0.3, 0.3])
ylabel({'Rate of change in'; 'Shoulder XY Rotation with Speed'; '(deg per (m/s))'})
% Statistics
[~,p_wlk_xy] = ttest(wlk_slopes_xy);
[~,p_run_xy] = ttest(run_slopes_xy);
text(0.75,-10,['p = ', num2str(p_wlk_xy)])
text(1.75,-10,['p = ', num2str(p_run_xy)])

subplot(2,2,3)
plot(sub_speeds{1,1}(1,1),Rotation{1,1}(1,2),'o','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',colors(1,:))
hold on;
for i = 1:size(subjects,1)
    for j = 1:size(Rotation{i,1},1)
        if dutyfactor{i,1}(j,1) >= 0.5 % walking
            plot(sub_speeds{i,1}(j,1),Rotation{i,1}(j,2),'o','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',colors(i,:))
        else % running
            plot(sub_speeds{i,1}(j,1),Rotation{i,1}(j,2),'s','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',colors(i,:))
        end
    end
    for k = 1:2 % splitting walking vs running
        if k == 1 % walking
            index = dutyfactor{i,1} >= 0.5;
        else % running
            index = dutyfactor{i,1} < 0.5;
        end
        % Linear Regresssion (assuming intercept should be Rot = 0)
        B_SH(:,i,k) = [ones(size(sub_speeds{i,1}(index,1))), sub_speeds{i,1}(index,1)]\Rotation{i,1}(index,2); % from Y = X[B0; B1], here B0 is the y-intercept and B1 is the slope
        y_reg_SH{i,k} = B_SH(2,i,k)*sub_speeds{i,1}(index,1) + B_SH(1,i,k); % y-values for the linear regression line
        plot(sub_speeds{i,1}(index,1),y_reg_SH{i,k},'-','Color',colors(i,:))
    end
end
ylabel({'Shoulder Rotation';'relative to hips (deg)'})
xlabel('Speed (m/s)')

% Determine average and standard deviation for shoulder rotation during 
% walking and running relative to hips 
wlk_slopes_SH = B_SH(2,:,1);
run_slopes_SH = B_SH(2,:,2);
Y_slope_SH_avg = [mean(wlk_slopes_SH,'omitnan'), mean(run_slopes_SH,'omitnan')];
subplot(2,2,4)
bar(X,Y_slope_SH_avg)
hold on;
er = errorbar(X,Y_slope_SH_avg,[std(wlk_slopes_SH,'omitnan'),std(run_slopes_SH,'omitnan')]);
er.Color = [0,0,0]; er.LineStyle = 'none';
X_rand_wlk = 0.1*(rand(size(wlk_slopes_SH))-0.5) + 1;
X_rand_run = 0.1*(rand(size(run_slopes_SH))-0.5) + 2;
plot(X_rand_wlk,wlk_slopes_SH,'.','Color',[1, 0.3, 0.3])
plot(X_rand_run,run_slopes_SH,'.','Color',[1, 0.3, 0.3])
ylabel({'Rate of change in'; 'Shoulder Rotation relative to Hips with Speed'; '(deg per (m/s))'})
% Statistics
[~,p_wlk_SH] = ttest(wlk_slopes_SH);
[~,p_run_SH] = ttest(run_slopes_SH);
text(0.75,23,['p = ', num2str(p_wlk_SH)])
text(1.75,23,['p = ', num2str(p_run_SH)])