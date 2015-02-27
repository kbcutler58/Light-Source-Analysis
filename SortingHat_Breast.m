%% Sorting Hat Script with Analysis
% Purpose: Load in Data for Light Source Wavelength and Fitting Analysis 
% Current Sorting Options
    % Data tag 1 = left breast
    % Data tag 2 = right breast
    % Data tag 3 = X Coordinate
    % Data tag 4 = Y Coordiante
    
    
% Load in Data
% clear
% clc
% cd('C:\Users\Kyle\Downloads\DOSI data for wv LG testing')
% 
% datafiles = dir('Breast Results\*waves01.mat');
% for i = 1:29
%     cd('C:\Users\Kyle\Downloads\DOSI data for wv LG testing')
%     dir_string = strcat('Breast Results\*waves',num2str(i,'%02d'),'.mat');
%     datafiles = dir(dir_string);
%     cd('Breast Results')
%     for totalcases = 1:length(datafiles)
%         load(datafiles(totalcases).name)
%         fieldlabel = strcat('Patient',num2str(totalcases,'%02d'));%,datafiles(totalcases+2).name(1:7));%end-4);
%         fieldlabel2 = strcat('waves',num2str(i,'%02d'));
%         test.(fieldlabel2).(fieldlabel)= PatientData;
%     end
% end
% 
%  disp('Data Loaded')
% 

%%

% for patientbatch = 1:26
%     patientSelect = patientbatch;

patientSelect = 1;


% Select for wavelengths
% waveSelect = 1:29;
waveSelect = 1;

waveSelection = zeros(1,29);
% waveSelection([1:5 7:9]) = 1; %1:29 %1 is original dataset
waveSelection(waveSelect) = 1; %1:29 %1 is original dataset
waveSelection = logical(waveSelection);
% Select for certain patients

patientSelection = zeros(1,24);
patientSelection(patientSelect) = 1; % 1:24
patientSelection = logical(patientSelection);

% Select for Points on Breast Patients
testcase1 = 1;
testcase2 = 0;

% % Exlusion Crit
excludePoints = 0;
if excludePoints == 1
    testcase3 = -40; testcase4 = 40; 
    testcase5 = -40; testcase6 = 40;
end

% Inclusion Crit
includePoints = 1;
if includePoints == 1
    testcase3 = 60; testcase4 = 60; 
    testcase5 = 70; testcase6 = 70;
%     testcase3 = -100; testcase4 = 100; 
%     testcase5 = -90; testcase6 = 0;
end

% includePoints = 1;
% if includePoints == 1
%     testcase3 = -60; testcase4 = -60; 
%     testcase5 = 40; testcase6 = 40;
% end

% Select for tissue type
% Not ready yet, possible to change test structure for correct tissue
TissueSelection = 1; % Breast = 1 Abs = 2 Calf = 3 Bicep = 4

fields = fieldnames(test);
selectedFields = fields(waveSelection);
clear selectedData
for i = 1:length(selectedFields)
    selectedData(i) = test.(selectedFields{i});
end

fields2 = fieldnames(selectedData);
selectedFields2 = fields2(patientSelection);
clear selectedData2
for i = 1:length(selectedFields2)
    for j = 1:length(selectedFields)
   selectedData2.(selectedFields{j}).(selectedFields2{i}) = selectedData(j).(selectedFields2{i});
%      selectedData2(j,i) = selectedData(j).(selectedFields2{i}); %Doesn't
%      work
    end
end

fields3 = fieldnames(selectedData2);
fields4 = fieldnames(selectedData2.(selectedFields{1}));
clear selectedDataFinal chromNew chromOld chromDif chromPerc muaOld muaNew muaOldFit wvs tags testwvs
muaOld = zeros(1,697);
muaNew = zeros(1,697);
wvs = zeros(1,697);
l = 1;
% Test the points to match point selection
for i = 1:length(selectedFields2)
    for j = 1:length(selectedFields)
        fields5 = fieldnames(selectedData2.(selectedFields{j}).(selectedFields2{i})); %Get fields for each patient
        for k = 1:length(fields5) % 1 : length of points
            if includePoints == 1
                if (selectedData2.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).tag(1) == testcase1) ...
                        && (selectedData2.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).tag(2) == testcase2) ...
                        && ((selectedData2.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).tag(3) >= testcase3) && (selectedData2.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).tag(3) <= testcase4))...
                        && ((selectedData2.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).tag(4) >= testcase5) && (selectedData2.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).tag(4) <= testcase6))
                    selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}) = selectedData2.(selectedFields{j}).(selectedFields2{i}).(fields5{k});
                    chromNew(l,:) =  selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).chromo(1,:);
                    chromOld(l,:) =  selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).chromo(2,:);
                    chromDif(l,:) =  selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).chromo(3,:);
                    chromPerc(l,:) =  selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).chromo(4,:);
                    muaSize = size(selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).muaSpectrum(2,:),2);
                    muaNew(l,1:muaSize) = selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).muaSpectrum(2,:);
                    muaOld(l,1:muaSize) = selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).muaSpectrum(3,:);
                    
                    if size(selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).muaSpectrum,1) > 3
                        muaOldFit(l,1:muaSize) = selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).muaSpectrum(4,:);
                    end
                    wvs(1,1:muaSize) = selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).muaSpectrum(1,:);
                    testwvs(l,:) = selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).wavelengths;
                    tags(l,:) = selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).tag;
                    l = l+1;
                end
            elseif excludePoints == 1
                if (selectedData2.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).tag(1) == testcase1) ...
                        && (selectedData2.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).tag(2) == testcase2) ...
                        && (((selectedData2.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).tag(3) <= testcase3) || (selectedData2.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).tag(3) >= testcase4))...
                        || ((selectedData2.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).tag(4) <= testcase5) || (selectedData2.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).tag(4) >= testcase6)))
                    selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}) = selectedData2.(selectedFields{j}).(selectedFields2{i}).(fields5{k});
                    chromNew(l,:) =  selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).chromo(1,:);
                    chromOld(l,:) =  selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).chromo(2,:);
                    chromDif(l,:) =  selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).chromo(3,:);
                    chromPerc(l,:) =  selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).chromo(4,:);
                    muaSize = size(selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).muaSpectrum(2,:),2);
                    muaNew(l,1:muaSize) = selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).muaSpectrum(2,:);
                    muaOld(l,1:muaSize) = selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).muaSpectrum(3,:); 
                    if size(selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).muaSpectrum,1) > 3
                        muaOldFit(l,1:muaSize) = selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).muaSpectrum(4,:);
                    end
                    wvs(1,1:muaSize) = selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).muaSpectrum(1,:);
                    testwvs(l,:) = selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).wavelengths;
                    tags(l,:) = selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).tag;
                    l = l+1;
                end
            end
        end
    end
end

results = ~exist('selectedDataFinal','var');
if results == 1
    disp('No Results Found')
end

% figure;
% scatter(tags(:,3),tags(:,4),500,chromDif(:,5),'fill'),colorbar %Hb02 Map
% figure;
% scatter(tags(:,3),tags(:,4),500,chromNew(:,5),'fill'),colorbar %Hb02 Map
% 

% wvs = selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).muaSpectrum(1,:);
% %% 
% range = 1:45;
% domain = find(muaOld(1,:) ~= 0);
% figure;
% titlestring = strcat('Patient ',num2str(patientSelect),' Wavelength Case ',num2str(waveSelect),' Broadband vs Discrete Point Fit');
% for r= 1:length(range);
%     hold on
%     plot(wvs,muaOld(r,:),'r',wvs,muaNew(r,:),'b'),legend('Broadband Fit','Discrete Point Fit'), title(titlestring),...
%         xlabel('wavelengths in nm'), ylabel('Mua in 1/mm');
% end

%% Plotting Mua Spectrum
% % Use multiple points on a patient

% % range = 1:45;
% range = 1:size(muaOld,1);
% % domain = find(muaOld(1,:) ~= 0)
% figure('units','normalized','outerposition',[0 0 1 1]);
% titlestring = strcat('Patient ',num2str(patientSelect),' Wavelength Case ',num2str(waveSelect),' Broadband vs Discrete Point Fit');
% for r= 1:length(range);
%     hold on
%     domain = find(muaOld(r,:) ~= 0);
%     plot(wvs(domain),muaOld(r,domain),'r',wvs(domain),muaNew(r,domain),'b'),legend('Broadband Fit','Discrete Point Fit'), title(titlestring),...
%         xlabel('wavelengths in nm'), ylabel('Mua in 1/mm');
%     [a,b,c] = intersect(testwvs(1,:),wvs);
%     redc = [1 1 1 1 1 1 1];
%     bluec = [0 0 0 0 0 0 0];
%     scatter(testwvs(1,:),muaOld(r,c),30,redc)
% %     scatter(testwvs(1,:),muaOld(r,c),30,redc,'fill')
%     scatter(testwvs(1,:),muaNew(r,c),30,bluec)
% %     scatter(testwvs(1,:),muaNew(r,c),30,bluec,'fill')
% end

% range = 1:45;
range = 1:size(muaOld,1);
% domain = find(muaOld(1,:) ~= 0)
% figure('units','normalized','outerposition',[0 0 1 1]);
titlestring = strcat('Patient ',num2str(patientSelect),' Wavelength Case ',num2str(waveSelect),' Broadband vs Discrete Point Fit');
for r= 1:length(range);
    hold on
    domain = find(muaOld(r,:) ~= 0);
    plot(wvs(domain),muaNew(r,domain),'b',wvs(domain),muaOld(r,domain),'r',wvs(domain),muaOldFit(r,domain),'g'),title(titlestring),...
        xlabel('wavelengths in nm'), ylabel('Mua in 1/mm');
    leg1 = legend('Discrete Point Fit','Broadband Fit','Broadband Conc');
    set(leg1,'Location','Best')
    [a,b,c] = intersect(testwvs(1,:),wvs);
    redc = [1 1 1 1 1 1 1];
    bluec = [0 0 0 0 0 0 0];
    scatter(testwvs(1,:),muaOld(r,c),30,redc)
%     scatter(testwvs(1,:),muaOld(r,c),30,redc,'fill')
    scatter(testwvs(1,:),muaNew(r,c),30,bluec)
%     scatter(testwvs(1,:),muaNew(r,c),30,bluec,'fill')
    scatter(testwvs(1,:),muaOldFit(r,c),30,bluec+.5)
%     scatter(testwvs(1,:),muaOld(r,c),30,redc,'fill')
end

% cd('C:\Users\Kyle\Downloads\DOSI data for wv LG testing\Analysis Pictures')
% h = gcf;
% savestring = strcat('Mua Spec ',num2str(patientSelect),'case',num2str(waveSelect(1)),'_',num2str(waveSelect(end)),' R','.jpg');
% saveas(h,savestring);

%% Plotting Concentration Differences
% Use multiple wavelengths and single patient
%  
% totalPoints = size(fieldnames(selectedDataFinal.(fields3{1}).(fields4{1})),1);
%  
% titlestring = strcat('Patient ',num2str(patientSelect),' Case ',num2str(waveSelect(1)),':',num2str(waveSelect(end)),' Conc Changes');
% figure('units','normalized','outerposition',[0 0 1 1]);
% subplot(2,2,1)
% for i = 1:length(waveSelect)
%     hold on
%     wvcombo = repmat(i,totalPoints,1 );
%     scatter(wvcombo,chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),1))
%     title(strcat(titlestring,' HbO2'))
%     xlabel('Wavelength Test Case')
%     ylabel('HbO2 difference in microMoles')
% end
% subplot(2,2,2)
% for i = 1:length(waveSelect)
%     hold on
%     wvcombo = repmat(i,totalPoints,1 );
%     scatter(wvcombo,chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),2))
%     title(strcat(titlestring,' Hb'))
%     xlabel('Wavelength Test Case')
%     ylabel('Hb difference in microMoles')
% 
% end
% subplot(2,2,3)
% for i = 1:length(waveSelect)
%     hold on
%     wvcombo = repmat(i,totalPoints,1 );
%     scatter(wvcombo,chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),3))
%     title(strcat(titlestring,' Water Frac'))
%     xlabel('Wavelength Test Case')
%     ylabel('Water Frac difference')
% 
% end
% subplot(2,2,4)
% for i = 1:length(waveSelect)
%     hold on
%     wvcombo = repmat(i,totalPoints,1 );
%     scatter(wvcombo,chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),4))
%     title(strcat(titlestring,' Lipid Frac'))
%     xlabel('Wavelength Test Case')
%     ylabel('Lipid Frac difference')
% end

% cd('C:\Users\Kyle\Downloads\DOSI data for wv LG testing\Analysis Pictures')
% h = gcf;
% savestring = strcat('NoBase Conc Changes ',num2str(patientSelect),' R','.jpg');
% saveas(h,savestring);


%% Plotting Concentration Differences in %
% Use multiple wavelengths and single patient
% 
% titlestring = strcat('Patient ',num2str(patientSelect),' Case ',num2str(waveSelect(1)),':',num2str(waveSelect(end)),' Changes in %');
% figure('units','normalized','outerposition',[0 0 1 1]);
% subplot(2,2,1)
% for i = 1:length(waveSelect)
%     hold on
%     wvcombo = repmat(i,totalPoints,1 );
%     scatter(wvcombo,chromPerc(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),1))
%     title(strcat(titlestring,' HbO2'))
%     xlabel('Wavelength Test Case')
%     ylabel('HbO2 difference in %')
% end
% subplot(2,2,2)
% for i = 1:length(waveSelect)
%     hold on
%     wvcombo = repmat(i,totalPoints,1 );
%     scatter(wvcombo,chromPerc(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),2))
%     title(strcat(titlestring,' Hb'))
%     xlabel('Wavelength Test Case')
%     ylabel('Hb difference in %')
% 
% end
% subplot(2,2,3)
% for i = 1:length(waveSelect)
%     hold on
%     wvcombo = repmat(i,totalPoints,1 );
%     scatter(wvcombo,chromPerc(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),3))
%     title(strcat(titlestring,' Water Frac'))
%     xlabel('Wavelength Test Case')
%     ylabel('Water Frac difference in %')
% 
% end
% subplot(2,2,4)
% for i = 1:length(waveSelect)
%     hold on
%     wvcombo = repmat(i,totalPoints,1 );
%     scatter(wvcombo,chromPerc(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),4))
%     title(strcat(titlestring,' Lipid Frac'))
%     xlabel('Wavelength Test Case')
%     ylabel('Lipid Frac difference in %')
% end

% 
% cd('C:\Users\Kyle\Downloads\DOSI data for wv LG testing\Analysis Pictures')
% h = gcf;
% savestring = strcat('NoBase Conc Changes in Perc ',num2str(patientSelect),' R','.jpg');
% saveas(h,savestring);
% close all




% end