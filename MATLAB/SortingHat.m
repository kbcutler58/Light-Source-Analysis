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
PremenSelect = [9 10 12 13 14 19 22 24 25]; %Patient Selection
PremenTest1 = [0 0 1 1 0 1 1 0 0 0]; %Left Breast Selection
PremenTest2 = [1 1 0 0 1 0 0 1 1 1]; %Right Breast Selection

% % Test Criteria for Healthy Postmenopausal Breast
PostmenSelect = [5 6 8 11 15 16 17 18 20 21]; %Patient Selection
PostmenTest1 = [0 0 1 1 0 0 1 1 0 1]; %Left Breast Selection
PostmenTest2 = [1 1 0 0 1 1 0 0 1 0]; %Right Breast Selection

% % Test Criteria for Breast Tumor Points
TumorSelect = [1 2 3 4 7 14 16 17 18 23]; %Patient Selection
TumorTest1 = [1 1 0 0 1 1 1 0 0 1]; %Left Breast Selection
TumorTest2 = [0 0 1 1 0 0 0 1 1 0]; %Right Breast Selection
TumorTest3 = [90 40 -80 -60 -20 -20 50 -80 -80 -100]; %X Lower Coordinate
TumorTest4 = [90 40 -80 -40 -10 -10 60 -70 -70 -80]; %X Upper Bound
TumorTest5 = [50 -40 50 30 60 80 -10 30 50 40 ]; %Y Lower Bound
TumorTest6 = [70 -20 70 30 70 100 0 40 50 40 ]; %Y Upper Bound

% Needed for increments of nonsequential patients
testcase_index = 0; %Select index of test criteria

% % Loop for sequential wavelength tests
% for waveSelect = 1:29;

% % Patient Selction
for patientbatch = 2; % Select Specific Patients% for patientbatch = [1 2 3 6 7 8 9 10 16] %Biceps Select Patients

PretrainSelect = [1 2 5 6 9 10 13 14 17 18 21 22 25 26];
PosttrainSelect = [3 4 7 8 11 12 15 16 19 20 23 24 27 28]; 

% for patientbatch = PremenSelect;
% for patientbatch = 1:26;
% for patientbatch = 1:21;
% for patientbatch = PretrainSelect;
% for patientbatch = PosttrainSelect;
% for patientbatch = 1:28;

% for patientbatch = TumorSelect
    patientSelect = patientbatch;

% for patientbatch = PostmenSelect % Loop for multipe patient cases
% patientSelect = PremenSelect(patientbatch) % Selection based on test criteria


% testcase_index = 10; % Select test criteria with specific index
testcase_index = testcase_index+1;

% % Patient Selection Groups
% patientSelect = 1:21; % Abdomen Patients
% patientSelect = 1:26; % Breast Patients
% patientSelect = 1:28; % Calf Patients
% patientSelect = 2; % Select individual patient

% % Select for wavelengths
% waveSelect = 1:29; % Use all wavelength cases
% waveSelect = waveSelect; % Use looped selection
waveSelect = 1:29;
% waveSelect = 1:17;

waveSelection = zeros(1,29);
% waveSelection([1:5 7:9]) = 1; %1:29 %1 is original dataset
waveSelection(waveSelect) = 1; %1:29 %1 is original dataset
waveSelection = logical(waveSelection);
% Select for certain patients

patientSelection = zeros(1,28); %Make logical array of highest number of patients
patientSelection(patientSelect) = 1; % 1:26
patientSelection = logical(patientSelection);

% Select for Points on Breast Patients
testcase1 = 0; % Left Breast
testcase2 = 1; % Right Breast

% Use Criteria for test case
% testcase1 = TumorTest1(testcase_index);
% testcase2 = TumorTest2(testcase_index);
% testcase1 = PostmenTest1(testcase_index);
% testcase2 = PostmenTest2(testcase_index);
% testcase1 = PremenTest1(testcase_index);
% testcase2 = PremenTest2(testcase_index);

% % Exlusion Crit (mostly for excluding areola on Breast)
excludePoints =0;
if excludePoints == 1
    testcase3 = -40; testcase4 = 40; 
    testcase5 = -40; testcase6 = 40;
end



% Inclusion Crit
includePoints = 1; % 1 to include a test range
if includePoints == 1
    testcase3 = -150; testcase4 = 150; 
    testcase5 = -150; testcase6 = 150;
%     testcase3 = -60; testcase4 = -60;
%     testcase5 = 40; testcase6 = 40;
end

% Specific Point Crit
specificPoint = 0; % Used for bounded point test cases
if specificPoint == 1
    testcase3 = TumorTest3(testcase_index);
    testcase4 = TumorTest4(testcase_index);
    testcase5 = TumorTest5(testcase_index);
    testcase6 = TumorTest6(testcase_index);
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
current_tag = [0; 0; 0; 0];
% Test the points to match point selection
for i = 1:length(selectedFields2)
    for j = 1:length(selectedFields)
        fields5 = fieldnames(selectedData2.(selectedFields{j}).(selectedFields2{i})); %Get fields for each patient
        for k = 1:length(fields5) % 1 : length of points
            if includePoints == 1 % Tag 1 & Tag 2 & (Tag3 & Tag4)
                if ((selectedData2.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).tag(1) == testcase1) ...
                        || (selectedData2.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).tag(2) == testcase2)) % ...
%                         && ((selectedData2.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).tag(3) >= testcase3) && (selectedData2.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).tag(3) <= testcase4))...
%                         && ((selectedData2.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).tag(4) >= testcase5) && (selectedData2.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).tag(4) <= testcase6))%% && (sum((selectedData2.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).tag(1:4) == current_tag)) < 4 )
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
%                     current_tag = tags(l,1:4)';
                    l = l+1;
                end
            elseif excludePoints == 1 % Tag 1 | Tag 2 & (Tag3 | Tag4)
                if ((selectedData2.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).tag(1) == testcase1) ...
                        || (selectedData2.(selectedFields{j}).(selectedFields2{i}).(fields5{k}).tag(2) == testcase2)) ...
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
%     for range = 1:totalPoints
figure;
       muaspec_plot 
%     end        
end

% close all
% muaspec_plot % Plot absorption spectra

end

% Appendix_Charter
% sum(Point_run)
% end

% % Calculate Chromophore Concentration Matrix
% clear TotalChrom
% TotalChrom(1,1) = sum((HbO2_run(1,:).*Point_run),2)./sum(Point_run);
% TotalChrom(2,1) = sum((HbO2_run(4,:).*Point_run),2)./sum(Point_run);
% TotalChrom(3,1) = min(HbO2_run(2,:));
% TotalChrom(4,1) = max(HbO2_run(3,:));
% 
% TotalChrom(1,2) = sum((Hb_run(1,:).*Point_run),2)./sum(Point_run);
% TotalChrom(2,2) = sum((Hb_run(4,:).*Point_run),2)./sum(Point_run);
% TotalChrom(3,2) = min(Hb_run(2,:));
% TotalChrom(4,2) = max(Hb_run(3,:));
% 
% TotalChrom(1,3) = sum((Water_run(1,:).*Point_run),2)./sum(Point_run);
% TotalChrom(2,3) = sum((Water_run(4,:).*Point_run),2)./sum(Point_run);
% TotalChrom(3,3) = min(Water_run(2,:));
% TotalChrom(4,3) = max(Water_run(3,:));
% 
% TotalChrom(1,4) = sum((Lipid_run(1,:).*Point_run),2)./sum(Point_run);
% TotalChrom(2,4) = sum((Lipid_run(4,:).*Point_run),2)./sum(Point_run);
% TotalChrom(3,4) = min(Lipid_run(2,:));
% TotalChrom(4,4) = max(Lipid_run(3,:));

