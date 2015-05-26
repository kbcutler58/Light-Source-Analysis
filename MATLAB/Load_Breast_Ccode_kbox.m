%% Load in Calf Measurement Analysis
% Purpose: Load in Data for Light Source Wavelength and Fitting Analysis 
   
    
%% Load in Data
% clear
% clc
% cd('C:\Users\Kyle\Downloads\DOSI data for wv LG testing')
% cd('C:\Users\Kyle\Downloads\c chromophore results\c chromophore results');
cd('C:\Users\Kyle\Downloads\breast4wvnew');
% datafiles = dir('breast6wv/*waves01*.mat');
datafiles = dir('breast4wvnew/*waves01*.mat');
% datafiles = dir('breast4wv/*waves01*.mat');
for i = 1:1
%     cd('C:\Users\Kyle\Downloads\c chromophore results\c chromophore results');
    cd('C:\Users\Kyle\Downloads\breast4wvnew');
%     dir_string = strcat('breast6wv/*waves',num2str(i,'%02d'),'*.mat');
%     dir_string = strcat('breast4wv/*waves',num2str(i,'%02d'),'*.mat');
    dir_string = strcat('breast4wvnew/*waves',num2str(i,'%02d'),'*.mat');
    datafiles = dir(dir_string);
    cd('breast4wvnew')
%     cd('breast6wv')
%     cd('breast4wv')
    for totalcases = 1:length(datafiles)
        load(datafiles(totalcases).name)
        fieldlabel = strcat('Patient',num2str(totalcases,'%02d'));%,datafiles(totalcases+2).name(1:7));%end-4);
        fieldlabel2 = strcat('waves',num2str(i,'%02d'));
        test2.(fieldlabel2).(fieldlabel)= PatientData;
    end
end

 disp('Data Loaded')

