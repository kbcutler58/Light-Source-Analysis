%% Load in Breast Measurement Analysis
% Purpose: Load in Data for Light Source Wavelength and Fitting Analysis 

%% Load in Data
clear
clc
cd('C:\Users\Kyle\Downloads\DOSI data for wv LG testing')

datafiles = dir('Breast Results 4wv\*waves01.mat');
for i = 1:17
    cd('C:\Users\Kyle\Downloads\DOSI data for wv LG testing')
    dir_string = strcat('Breast Results 4wv\*waves',num2str(i,'%02d'),'.mat');
    datafiles = dir(dir_string);
    cd('Breast Results 4wv')
    for totalcases = 1:length(datafiles)
        load(datafiles(totalcases).name)
        fieldlabel = strcat('Patient',num2str(totalcases,'%02d'));%,datafiles(totalcases+2).name(1:7));%end-4);
        fieldlabel2 = strcat('waves',num2str(i,'%02d'));
        test.(fieldlabel2).(fieldlabel)= PatientData;
    end
end

 disp('Data Loaded')
