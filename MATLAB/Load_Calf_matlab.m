%% Load in Calf Measurement Analysis
% Purpose: Load in Data for Light Source Wavelength and Fitting Analysis 
   
    
%% Load in Data
clear
clc
cd('C:\Users\Kyle\Downloads\DOSI data for wv LG testing')

datafiles = dir('Calf Results/*waves01*.mat');
for i = 1:29
    cd('C:\Users\Kyle\Downloads\DOSI data for wv LG testing')
    dir_string = strcat('Calf Results/*waves',num2str(i,'%02d'),'*.mat');
    datafiles = dir(dir_string);
    cd('Calf Results')
    for totalcases = 1:length(datafiles)
        load(datafiles(totalcases).name)
        fieldlabel = strcat('Patient',num2str(totalcases,'%02d'));%,datafiles(totalcases+2).name(1:7));%end-4);
        fieldlabel2 = strcat('waves',num2str(i,'%02d'));
        test.(fieldlabel2).(fieldlabel)= PatientData;
    end
end

 disp('Data Loaded')

