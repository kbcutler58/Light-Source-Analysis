%% Sorting Hat Script with Analysis
% Sorting Data for Light Source Wavelength and Fitting Analysis 
% Current Sorting Options
    % Data tag 1 = left 
    % Data tag 2 = right
    % Data tag 3 = X Coordinate
    % Data tag 4 = Y Coordiante

%% Declare test criteria, declare variables

% Running averages for multiple patient test cases 
HbO2_run = [];
Hb_run = [];
Water_run = [];
Lipid_run = [];
Point_run = [];

% Test Criteria for Healthy Premenopausal Breast
% PremenSelect = [9 10 12 13 14 19 22 24 25]; %Patient Selection
% PremenTest1 = [0 0 1 1 0 1 1 0 0 0]; %Left Breast Selection
% PremenTest2 = [1 1 0 0 1 0 0 1 1 1]; %Right Breast Selection

% % Test Criteria for Healthy Postmenopausal Breast
% PostmenSelect = [5 6 8 11 15 16 17 18 20 21]; %Patient Selection
% PostmenTest1 = [0 0 1 1 0 0 1 1 0 1]; %Left Breast Selection
% PostmenTest2 = [1 1 0 0 1 1 0 0 1 0]; %Right Breast Selection

% % Test Criteria for Breast Tumor Points
% TumorSelect = [1 2 3 4 7 14 16 17 18 23]; %Patient Selection
% TumorTest1 = [1 1 0 0 1 1 1 0 0 1]; %Left Breast Selection
% TumorTest2 = [0 0 1 1 0 0 0 1 1 0]; %Right Breast Selection
% TumorTest3 = [90 40 -80 -60 -20 -20 50 -80 -80 40]; X Lower Coordinate
% TumorTest4 = [90 40 -80 -40 -10 -10 60 -70 -70 40]; X Upper Bound
% TumorTest5 = [50 -40 50 30 60 80 -10 30 50 -100]; Y Lower Bound
% TumorTest6 = [70 -20 70 30 70 100 0 40 50 -80]; Y Upper Bound

% Needed for increments of nonsequential patients
% testcase_index = 0; %Select index of test criteria

% % Loop for sequential wavelength tests
% for waveSelect = 1:29;

% % Patient Selction
% for patientbatch = [16 19]; % Select Specific Patients
% for patientbatch = [1 2 3 6 7 8 9 10 16] %Biceps Select Patients

% for patientbatch = PostmenSelect % Loop for multipe patient cases
% patientSelect = PremenSelect(patientbatch) % Selection based on test criteria

% testcase_index = 3 % Select test criteria with specific index
% testcase_index = testcase_index+1;

% % Patient Selection Groups
% patientSelect = 1:21; % Abdomen Patients
% patientSelect = 1:26; % Breast Patients
% patientSelect = 1:28; % Calf Patients
% patientSelect = 1; % Select individual patient

% % Select for wavelengths
% waveSelect = 1:29; % Use all wavelength cases
% waveSelect = waveSelect; % Use looped selection
waveSelect = 1;

waveSelection = zeros(1,29);
% waveSelection([1:5 7:9]) = 1; %1:29 %1 is original dataset
waveSelection(waveSelect) = 1; %1:29 %1 is original dataset
waveSelection = logical(waveSelection);
% Select for certain patients

patientSelection = zeros(1,26); %Make logical array of highest number of patients
patientSelection(patientSelect) = 1; % 1:26
patientSelection = logical(patientSelection);

% Select for Points on Breast Patients
% testcase1 = 0; % Left Breast
% testcase2 = 0; % Right Breast

% Use Criteria for test case
% testcase1 = TumorTest1(testcase_index);
% testcase2 = TumorTest2(testcase_index);
% testcase1 = PostmenTest1(testcase_index);
% testcase2 = PostmenTest2(testcase_index);
% testcase1 = PremenTest1(testcase_index);
% testcase2 = PremenTest2(testcase_index);

% % Exlusion Crit (mostly for excluding areola on Breast)
excludePoints =1;
if excludePoints == 1
    testcase3 = -30; testcase4 = 30; 
    testcase5 = -30; testcase6 = 30;
end

% Specific Point Crit
specificPoint = 0; % Used for bounded point test cases
if specificPoint == 1
    testcase3 = TumorTest3(testcase_index);
    testcase4 = TumorTest4(testcase_index);
    testcase5 = TumorTest5(testcase_index);
    testcase6 = TumorTest6(testcase_index);
end

% Inclusion Crit
includePoints = 0; % 1 to include a test range
if includePoints == 1
    testcase3 = -150; testcase4 = 150; 
    testcase5 = -150; testcase6 = 150;
%     testcase3 = -60; testcase4 = -60;
%     testcase5 = 40; testcase6 = 40;
end

% Select for tissue type
% Not ready yet, possible to change test structure for correct tissue
TissueSelection = 1; % Breast = 1 Abs = 2 Calf = 3 Bicep = 4

% Select for wavelengths and load into data structure
fields = fieldnames(test);
selectedFields = fields(waveSelection);
clear selectedData
for i = 1:length(selectedFields)
    selectedData(i) = test.(selectedFields{i});
end

% Select for patients in test criteria
fields2 = fieldnames(selectedData);
selectedFields2 = fields2(patientSelection);
clear selectedData2
for i = 1:length(selectedFields2)
    for j = 1:length(selectedFields)
   selectedData2.(selectedFields{j}).(selectedFields2{i}) = selectedData(j).(selectedFields2{i});
    end
end

% Select for points by using data tags / inclusion or exclusion
fields3 = fieldnames(selectedData2);
fields4 = fieldnames(selectedData2.(selectedFields{1}));
clear selectedDataFinal chromNew chromOld chromDif chromPerc muaOld muaNew muaOldFit wvs tags testwvs
muaOld = zeros(1,697); % Used for different spec measurement wvs
muaNew = zeros(1,697); % Used for different spec measurement wvs
muaOldFit = zeros(1,697); % Used for different spec measurement wvs
wvs = zeros(1,697); % Used for different spec measurement wvs

l = 1; % Variable used to extract all points in single matrix

% Test the points to match point selection
for i = 1:length(selectedFields2)
    for j = 1:length(selectedFields)
        fields5 = fieldnames(selectedData2.(selectedFields{j}).(selectedFields2{i})); %Get fields for each patient
        for k = 1:length(fields5) % 1 : length of points
            if includePoints == 1 % Tag 1 & Tag 2 & (Tag3 & Tag4)
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
%                     muaOldFit(l,1:muaSize) = selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).muaSpectrum(4,:);
                    wvs(1,1:muaSize) = selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).muaSpectrum(1,:);
                    testwvs(l,:) = selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).wavelengths;
                    tags(l,:) = selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).tag;
                    l = l+1;
                end
            elseif excludePoints == 1 % Tag 1 & Tag 2 & (Tag3 || Tag4)
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
%                     muaOldFit(l,1:muaSize) = selectedDataFinal.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).muaSpectrum(4,:);
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
else    
    chromo_plot % Plot conc vs wavelengths
    HbO2_run = [HbO2_run HbO2];
    Hb_run = [Hb_run Hb];
    Water_run = [Water_run Water];
    Lipid_run = [Lipid_run Lipid];
    Point_run = [Point_run totalPoints];
end

% close all
% muaspec_plot % Plot absorption spectra

% end
% end

% Calculate Chromophore Concentration Matrix
clear TotalChrom
TotalChrom(1,1) = sum((HbO2_run(1,:).*Point_run),2)./sum(Point_run);
TotalChrom(2,1) = sum((HbO2_run(4,:).*Point_run),2)./sum(Point_run);
TotalChrom(3,1) = min(HbO2_run(2,:));
TotalChrom(4,1) = max(HbO2_run(3,:));

TotalChrom(1,2) = sum((Hb_run(1,:).*Point_run),2)./sum(Point_run);
TotalChrom(2,2) = sum((Hb_run(4,:).*Point_run),2)./sum(Point_run);
TotalChrom(3,2) = min(Hb_run(2,:));
TotalChrom(4,2) = max(Hb_run(3,:));

TotalChrom(1,3) = sum((Water_run(1,:).*Point_run),2)./sum(Point_run);
TotalChrom(2,3) = sum((Water_run(4,:).*Point_run),2)./sum(Point_run);
TotalChrom(3,3) = min(Water_run(2,:));
TotalChrom(4,3) = max(Water_run(3,:));

TotalChrom(1,4) = sum((Lipid_run(1,:).*Point_run),2)./sum(Point_run);
TotalChrom(2,4) = sum((Lipid_run(4,:).*Point_run),2)./sum(Point_run);
TotalChrom(3,4) = min(Lipid_run(2,:));
TotalChrom(4,4) = max(Lipid_run(3,:));

%{
%%
% end

% totalPoints = size(fieldnames(selectedDataFinal.(fields3{1}).(fields4{1})),1);
% disp(totalPoints)
 
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
% 
% range = 1;%
% % :45;
% % range = 1:size(muaOld,1);
% % domain = find(muaOld(1,:) ~= 0)
% % figure('units','normalized','outerposition',[0 0 1 1]);
% titlestring = strcat('Patient ',num2str(patientSelect),' Wavelength Case ',num2str(waveSelect(1)),':',num2str(waveSelect(end)),' Broadband vs Discrete Point Fit');
% for r= 1:length(range);
%     hold on
%     domain = find(muaOld(r,:) ~= 0);
%     plot(wvs(domain),muaNew(r,domain),'b',wvs(domain),muaOld(r,domain),'r',wvs(domain),muaOldFit(r,domain),'g','LineWidth',1.5),title(titlestring),...
%         xlabel('wavelengths in nm'), ylabel('{\mu}_a (1/mm)');
%     leg1 = legend('Discrete Point Fit','Broadband Fit','Broadband Conc');
% %     han = gtext('Example String')
%     ax = gca;
% 
% %     set (ax,'FontName','Symbol')
%     set(leg1,'Location','Best')
%     
%     [a,b,c] = intersect(testwvs(1,:),wvs);
%     redc = [1 1 1 1 1 1 1];
%     bluec = [0 0 0 0 0 0 0];
%     scatter(testwvs(1,:),muaOld(r,c),30,redc)
%     scatter(testwvs(1,:),muaNew(r,c),30,bluec)
%     scatter(testwvs(1,:),muaOldFit(r,c),30,bluec+.5)
%     
% %     scatter(testwvs(1,:),muaOld(r,c),30,redc,'fill')
% %     scatter(testwvs(1,:),muaNew(r,c),30,bluec,'fill')
% %     scatter(testwvs(1,:),muaOldFit(r,c),30,bluec+.5,'fill')
% 
% end

%% Conc

% 
% totalPoints = size(fieldnames(selectedDataFinal.(fields3{1}).(fields4{1})),1);
% disp(totalPoints);
%  clear Hb HbO2 Water Lipid
% titlestring = strcat('Patient ',num2str(patientSelect(1)),':',num2str(patientSelect(end)),' Case ',num2str(waveSelect(1)),':',num2str(waveSelect(end)),' Conc Changes');
% figure('units','normalized','outerposition',[0 0 1 1]);
% subplot(2,2,1)
% for i = 1:length(waveSelect)
%     hold on
%     wvcombo = repmat(i,totalPoints,1 );
%     scatter(wvcombo,chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),1))
%     title(strcat(titlestring,' HbO2'))
%     xlabel('Wavelength Test Case')
%     ylabel('HbO2 difference in {\mu}M')
%     ax = gca;
%     set(ax,'YLim',[(-max(abs(chromDif(:,1))))-1 max(abs(chromDif(:,1)))+1])
%     disp(strcat('Case ',num2str(i,'%02d')));
%     HbO2(1,i) = mean(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),1));
%     HbO2(2,i) = min(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),1));
%     HbO2(3,i) = max(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),1));
%     HbO2(4,i) = std(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),1));
% 
% %     hline.Color = 'black';
% %     set(ax,'YGrid','on')
% %     set(ax,'YTick',0)
% %     set(ax,'YMinorTick')
% %     set(ax,'Xtick',[])
% %     grid on
%     
% end
%     hline = refline([0 0]);
%     set(hline,'Color','k')
%     
% subplot(2,2,2);
% for i = 1:length(waveSelect)
%     hold on
%     wvcombo = repmat(i,totalPoints,1 );
%     scatter(wvcombo,chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),2))
%     title(strcat(titlestring,' Hb'))
%     xlabel('Wavelength Test Case')
%     ylabel('Hb difference in {\mu}M')
%     ax = gca;
%     set(ax,'YLim',[(-max(abs(chromDif(:,2))))-1 max(abs(chromDif(:,2)))+1])
%     Hb(1,i) = mean(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),2));
%     Hb(2,i) = min(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),2));
%     Hb(3,i) = max(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),2));
%     Hb(4,i) = std(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),2));
% %     set(ax,'YGrid','on')
% %     set(ax,'YTick',0)
% %     set(ax,'YMinorTick')
% %     set(ax,'Xtick',[])
% %     grid on
% 
% end
%     hline = refline([0 0]);
%     set(hline,'Color','k')
%     
% subplot(2,2,3)
% for i = 1:length(waveSelect)
%     hold on
%     wvcombo = repmat(i,totalPoints,1 );
%     scatter(wvcombo,chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),3))
%     title(strcat(titlestring,' Water Frac'))
%     xlabel('Wavelength Test Case')
%     ylabel('Water Frac difference')
%     ax = gca;
%     set(ax,'YLim',[(-max(abs(chromDif(:,3))))-1 max(abs(chromDif(:,3)))+1])
%     Water(1,i) = mean(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),3));
%     Water(2,i) = min(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),3));
%     Water(3,i) = max(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),3));
%     Water(4,i) = std(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),3));
% %     set(ax,'YGrid','on')
% %     set(ax,'YTick',0)
% %     set(ax,'YMinorTick')
% %     set(ax,'Xtick',[])
% %     grid on
% 
% end
% 
%     hline = refline([0 0]);
%     set(hline,'Color','k')
%     
% subplot(2,2,4)
% for i = 1:length(waveSelect)
%     hold on
%     wvcombo = repmat(i,totalPoints,1 );
%     scatter(wvcombo,chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),4))
%     title(strcat(titlestring,' Lipid Frac'))
%     xlabel('Wavelength Test Case')
%     ylabel('Lipid Frac difference')
%     ax = gca;
%     set(ax,'YLim',[(-max(abs(chromDif(:,4))))-1 max(abs(chromDif(:,4)))+1])
%     Lipid(1,i) = mean(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),4));
%     Lipid(2,i) = min(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),4));
%     Lipid(3,i) = max(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),4));
%     Lipid(4,i) = std(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),4));
%     
% %     set(ax,'YGrid','on')
% %     set(ax,'YTick',0)
% %     set(ax,'YMinorTick')
% %     set(ax,'Xtick',[])
% %     grid on
% end
%     hline = refline([0 0]);
%     set(hline,'Color','k')
%     
%% Plotting Concentration Differences
% Use multiple wavelengths and single patient
 
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
% 
% cd('C:\Users\Kyle\Downloads\DOSI data for wv LG testing\Analysis Pictures')
% h = gcf;
% savestring = strcat('Conc Changes ',num2str(patientSelect),' Abs','.jpg');
% saveas(h,savestring);


%% Plotting Concentration Differences in %
% Use multiple wavelengths and single patient

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
% savestring = strcat('Conc Changes in Perc ',num2str(patientSelect),' Abs','.jpg');
% saveas(h,savestring);
% close all




% end
%}