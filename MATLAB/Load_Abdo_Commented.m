%% Load in Abdomen Measurement Analysis
% Purpose: Load in Data for Light Source Wavelength and Fitting Analysis 
   
%% Load in Data
clear
clc
cd('C:\Users\Kyle\Downloads\DOSI data for wv LG testing') % Data Directory 

datafiles = dir('Ab Results 4wv/*waves01*.mat'); %Subfolder, results to get numer of files

for i = 1:17
    cd('C:\Users\Kyle\Downloads\DOSI data for wv LG testing') % Go into subfolder
    dir_string = strcat('Ab Results 4wv/*waves',num2str(i,'%02d'),'*.mat'); %Find files matching wavelength case
    datafiles = dir(dir_string); %file names
    cd('Ab Results 4wv') % go to specific folder
    for totalcases = 1:length(datafiles)
        load(datafiles(totalcases).name) % Load in each data file
        fieldlabel = strcat('Patient',num2str(totalcases,'%02d')); %Stucture tag 1: Patient %,datafiles(totalcases+2).name(1:7));%end-4);
        fieldlabel2 = strcat('waves',num2str(i,'%02d')); % Structure tag 2: Wavelength case
        test.(fieldlabel2).(fieldlabel)= PatientData; % Create structure using tags
    end
end

 disp('Data Loaded')
