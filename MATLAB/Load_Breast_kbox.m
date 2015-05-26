%% Load in Breast Measurement Analysis
% Purpose: Load in Data for Light Source Wavelength and Fitting Analysis 

%% Load in Data
% clear
% clc
cd('C:\Users\Kyle\Downloads\DOSI data for wv LG testing')

% datafiles = dir('Breast Results 4wv\*waves01.mat');
datafiles = dir('Breast Results 4wv new\*waves01.mat');
% datafiles = dir('Breast Results\*waves01.mat');
for i = 1:1
    cd('C:\Users\Kyle\Downloads\DOSI data for wv LG testing')
    dir_string = strcat('Breast Results 4wv new\*waves',num2str(i,'%02d'),'.mat');
%     dir_string = strcat('Breast Results 4wv\*waves',num2str(i,'%02d'),'.mat');
%         dir_string = strcat('Breast Results\*waves',num2str(i,'%02d'),'.mat');
    datafiles = dir(dir_string);
%     cd('Breast Results 4wv')
    cd('Breast Results 4wv new')
%      cd('Breast Results')
    for totalcases = 1:length(datafiles)
        load(datafiles(totalcases).name)
        fieldlabel = strcat('Patient',num2str(totalcases,'%02d'));
        fieldlabel2 = strcat('waves',num2str(i,'%02d'));
        test.(fieldlabel2).(fieldlabel)= PatientData;
    end
end

 disp('Data Loaded')
