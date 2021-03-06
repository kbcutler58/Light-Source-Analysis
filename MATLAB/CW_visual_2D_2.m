% CW 2D Visualization Program
% Kyle Cutler 7/31/2014
% Last Updated 9/11/2014

%{
Scripted design to process and display data taken from CW/Tracking System
%}
clear
clc

%% Load in Data

%txt file output path
pathName = 'C:\Users\Kyle\Downloads';
outputFile = 'test_txt.txt';

cd(pathName) % path folder
data = load(outputFile); % import data
%% Specify Data Format

label.time = 1;
label.x = 2;
label.y = 3;
label.yaw = 4;
label.pitch = 5;
label.roll = 6;
label.button = 7;
label.squal = 8;
% label.random = 9;
label.waveform1 = 10:84;
label.waveform2 = 85:159;

data_time = data(:,label.time);
data_x = data(:,label.x);
data_y = data(:,label.y);
data_yaw = deg2rad(data(:,label.yaw));
data_pitch = deg2rad(data(:,label.pitch));
data_roll = deg2rad(data(:,label.roll));
data_squal = data(:,label.squal);
% data_random = data(:,label.random);
data_button = data(:,label.button);
data_waveform1 = data(:,label.waveform1);
data_waveform2 = data(:,label.waveform2);

clear data outputFile pathName label
%% Find Peak Power in Waveforms

FFT_Resolution = 14;
peak_power = zeros(size(data_waveform1,1),4);

for i=1:size(data_waveform1,1)
   [peak_power(i,1),peak_power(i,2)]=CWFFT3(data_waveform1(i,:),FFT_Resolution);
   [peak_power(i,3),peak_power(i,4)]=CWFFT3(data_waveform2(i,:),FFT_Resolution);
end

% data = data3; % renaming probe output
% peaks = peaks3; % renaming fft output
% colordata = [peaks(:,1),peaks(:,3)]; % fft data to color output

%% Plot raw x,y path with color

% peak_power(:,1) = 20k peak
% peak_power(:,2) = 20k freq
% peak_power(:,3) = 11k peak
% peak_power(:,4) = 11k freq

% plot integrated x,y output with color from peak data

% %Plot 20khz and raw x and y
% scatter(cumsum(data_x),cumsum(data_y),150,-peak_power(:,3)); colorbar, axis square

% %Plot 11khz and raw x and y
% scatter(cumsum(data_x),cumsum(data_y),150,-peak_power(:,1)); colorbar, axis square

%% Correct for rotation changes of probe
% data_x = raw x data
% data_y = raw y data
% data_yaw = convert yaw to rads

data_yaw_offset = data_yaw - data_yaw(1); % offset angle to zero
x_location = 0; % initial x location
y_location = 0; % initial y location
x_location_offset = 0; % initial x with no offset
y_location_offset = 0; %initial y with no offset
corrected_path = zeros(size(data_x,1),2);
corrected_path_offset = zeros(size(data_x,1),2);

% No angle offset integration
for i=1:length(data_x)
    % Calculate rotation matrix
    correction_factor = [ cos(data_yaw(i)) -sin(data_yaw(i)); sin(data_yaw(i)) cos(data_yaw(i))];
    % Matrix for raw x and y
    path_matrix = [data_x(i),data_y(i)];
    % Multiply x,y matrix and rotation matrix
    path_temp = path_matrix*correction_factor;
    % Add temp path to previous location
    corrected_path(i,:) = [x_location + path_temp(1) , y_location + path_temp(2)];
    % Update Location of x and y
    x_location = corrected_path(i,1);
    y_location = corrected_path(i,2);
end

% Angle offset integration
for i=1:length(data_x)
    % Calculate rotation matrix 
    correction_factor = [ cos(data_yaw_offset(i)) -sin(data_yaw_offset(i)); sin(data_yaw_offset(i)) cos(data_yaw_offset(i))];
    % Matrix for raw x and y
    path_matrix = [data_x(i),data_y(i)];
    % Multiply x,y matrix and rotation matrix
    path_temp = path_matrix*correction_factor;
    % Add temp path to previous location
    corrected_path_offset(i,:) = [x_location_offset + path_temp(1) , y_location_offset + path_temp(2)];
    % Update Location of x and y 
    x_location_offset = corrected_path_offset(i,1);
    y_location_offset = corrected_path_offset(i,2);
end

%% Plotting 
markerSize = 150;
% axis([xmin xmax ymin ymax]
newaxis = ([-10000 10000 -10000 10000]);

% % Plot x y raw path
% close all
% scatter(cumsum(data_x),cumsum(data_y),150,-peak_power(:,3)); colorbar,
% % axis(newaxis); 
% axis equal
% title('2D Raw X and Y Output');

% Plot x y angle corrected path with offset (11khz peak)
figure;
scatter(corrected_path_offset(:,1),corrected_path_offset(:,2),markerSize,-peak_power(:,3)), colorbar
title('2D Path Angle Corrected Offset');
% axis(newaxis); 
axis equal

% % Plot x y angle corrected path with offset (20khz peak)
% figure;
% scatter(corrected_path_offset(:,1),corrected_path_offset(:,2),markerSize,-peak_power(:,3)), colorbar
% title('2D Path Angle Corrected Offset');
% % axis(newaxis); 
% axis equal

%% Plot overlay 
% 
% close all
% hold on
% title('2D Path Angle Corrected')
% xlabel('X Axis in pixels')
% ylabel('Y Axis in pixels')
% scatter(cumsum(data(:,2)),cumsum(data(:,3)),markerSize,-colordata(:,1)), colorbar;
% scatter(corrected_path2(:,1),corrected_path2(:,2),markerSize, -colordata(:,1)), colorbar;
% legend('original path','corrected path')
% % axis(newaxis); 
% axis equal
