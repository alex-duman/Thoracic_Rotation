% This script is intended to read in the treadmill data file you choose,
% check via plotting whether the data was zeroed and then prompt you to
% zero the data if it is not by selecting ranges over which the values
% should be zeroed.

% Select File
[file,path] = uigetfile('*Trd_Filt.csv','Select Treadmill file to check zeroing.');

% Open the file
Data = readtable([path file]);

% Instructions
m1 = msgbox({'In the following graph, select a point somewhere on the graph when ';
            'both the front and rear force plates appear to have brief periods ';
            'that should be zeroed (flat horizontal stretches). This will be used';
            'to generate a new zoomed-in graph around that region.';
            '';
            'Hit Return/Enter when done.'});
uiwait(m1)

fig1a = figure('Name','Pick point where both plates have zero condition');
subplot(2,1,1)
plot(Data.Front_Fx,'r-')
hold on;
plot(Data.Front_Fy,'g-')
plot(Data.Front_Fz,'b-')
ylabel('Force (N)')
legend('F_x','F_y','F_z','Location','east')

subplot(2,1,2)
plot(Data.Rear_Fx,'r-')
hold on;
plot(Data.Rear_Fy,'g-')
plot(Data.Rear_Fz,'b-')
ylabel('Force (N)')
xlabel('dimension of force vector')
[x,~] = ginput(1);
close(fig1a);
x = round(x);

if x < 2501 % if within first 2.5 seconds
    x = 2501; % recenter so we get first 5 second window
elseif x > length(Data.Front_Fx)-2500 % if within last 2.5 seconds
    x = length(Data.Front_Fx)-2500; % recenter to get last 5 second window
end

% Instructions
m2 = msgbox({'In the following graph, select the ranges along the Fz (vertical force)';
            'which correspond to when the value should be zero (subject is not on';
            'the force plate. Start from left-most side of x-axis and work toward';
            'the right, only need to sample a few timeframes when it should be zero.';
            '';
            'Hit Return/Enter when done.'});
uiwait(m2)

% Display Front Force Plate Data
fig1 = figure('Name','Front Force Plate');
plot(Data.time(x-2500:x+2500),Data.Front_Fx(x-2500:x+2500),'r-')
hold on;
plot(Data.time(x-2500:x+2500),Data.Front_Fy(x-2500:x+2500),'g-')
plot(Data.time(x-2500:x+2500),Data.Front_Fz(x-2500:x+2500),'b-')
title('5 second region of interest')
ylabel('Force (N)')
xlabel('time (s)')
legend('F_x','F_y','F_z','Location','east')

[X_front,~] = ginput; % collecting points of interest
close(fig1)

% Display Rear Force Plate Data
fig2 = figure('Name','Rear Force Plate');
plot(Data.time(x-2500:x+2500),Data.Rear_Fx(x-2500:x+2500),'r-')
hold on;
plot(Data.time(x-2500:x+2500),Data.Rear_Fy(x-2500:x+2500),'g-')
plot(Data.time(x-2500:x+2500),Data.Rear_Fz(x-2500:x+2500),'b-')
title('5 second region of interest')
ylabel('Force (N)')
xlabel('time (s)')
legend('F_x','F_y','F_z','Location','east')

[X_rear,~] = ginput;
close(fig2)

if length(X_front) > 1 % only need to rezero if points were selected
    % Determining Average Values for regions in Front Plate
    row = 1;
    for i = 1:2:length(X_front)
        Front_Fx_offsets(row,1) = mean(Data.Front_Fx(round(X_front(i)*1000):round(X_front(i+1)*1000)));
        Front_Fy_offsets(row,1) = mean(Data.Front_Fy(round(X_front(i)*1000):round(X_front(i+1)*1000)));
        Front_Fz_offsets(row,1) = mean(Data.Front_Fz(round(X_front(i)*1000):round(X_front(i+1)*1000)));
        row = row + 1;
    end
    Front_Fx_offset = mean(Front_Fx_offsets);
    Front_Fy_offset = mean(Front_Fy_offsets);
    Front_Fz_offset = mean(Front_Fz_offsets);
    % Rezeroing Values
    Data.Front_Fx = Data.Front_Fx - Front_Fx_offset;
    Data.Front_Fy = Data.Front_Fy - Front_Fy_offset;
    Data.Front_Fz = Data.Front_Fz - Front_Fz_offset;
end

if length(X_rear) > 1 % only need to rezero if points were selected
    % Determining Average Values for regions in Rear Plate
    row = 1;
    for i = 1:2:length(X_rear)
        Rear_Fx_offsets(row,1) = mean(Data.Rear_Fx(round(X_rear(i)*1000):round(X_rear(i+1)*1000)));
        Rear_Fy_offsets(row,1) = mean(Data.Rear_Fy(round(X_rear(i)*1000):round(X_rear(i+1)*1000)));
        Rear_Fz_offsets(row,1) = mean(Data.Rear_Fz(round(X_rear(i)*1000):round(X_rear(i+1)*1000)));
        row = row + 1;
    end
    Rear_Fx_offset = mean(Rear_Fx_offsets);
    Rear_Fy_offset = mean(Rear_Fy_offsets);
    Rear_Fz_offset = mean(Rear_Fz_offsets);
    % Rezeroing Values
    Data.Rear_Fx = Data.Rear_Fx - Rear_Fx_offset;
    Data.Rear_Fy = Data.Rear_Fy - Rear_Fy_offset;
    Data.Rear_Fz = Data.Rear_Fz - Rear_Fz_offset;
end

% Final Plot to Check
fig3 = figure('Name','Check Rezeroed Force Data Before Saving');
subplot(2,1,1)
plot(Data.time,Data.Front_Fx,'r-')
hold on;
plot(Data.time,Data.Front_Fy,'g-')
plot(Data.time,Data.Front_Fz,'b-')
legend('F_x','F_y','F_z','Location','east')
ylabel('Force (N)')
xlabel('time (s)')
title('Front Plate')

subplot(2,1,2)
plot(Data.time,Data.Rear_Fx,'r-')
hold on;
plot(Data.time,Data.Rear_Fy,'g-')
plot(Data.time,Data.Rear_Fz,'b-')
xlabel('time (s)')
title('Rear Plate')
uiwait(fig3) % wait for user to close figure before prompting about saving data

% Asking user to save data or not to original file name
quest = 'Do you want to overwrite the Treadmill data with this new force data?';
answer = questdlg(quest,'Saving Option','Overwrite','Keep Original','Keep Original');
% Handling Responses
switch answer
    case 'Overwrite'
        writetable(Data, [path file]);
        disp('Data was successfully updated!')
    case 'Keep Original'
        disp('Original file has been kept!')
end

