%% Load in Breast Data from C Data
clear
clc

cd('C:\Users\Kyle\Documents\MATLAB')
Load_Breast_Ccode_kbox;
cd('C:\Users\Kyle\Documents\MATLAB')
Load_Breast_kbox_nonneg;
tissue = 'Breast 4wv kbox nonneg'
% Load_Breast_Ccode
% cd('C:\Users\Kyle\Documents\MATLAB')
% Load_Breast
% tissue='Breast 6wv'
% tissue='Breast 4wv'
%  
% Load_Calf_Ccode
% Load_Calf
% tissue='Calf 4wv';
% tissue='Calf 6wv';


%%
cd('C:\Users\Kyle\Documents\MATLAB\Analysis Pictures')
fid = fopen(strcat(tissue,'muaDiff.txt'),'w');
fid2 = fopen(strcat(tissue,'muaPerdiff.txt'),'w');
fid3 = fopen(strcat(tissue,'muaStd.txt'),'w');

HbO2_run = [];
Hb_run = [];
Water_run = [];
Lipid_run = [];
Point_run = [];
% Mua_run = [];
% Mua_old = [];

for i = 1:26 % Breast
% for i = 1:28 % Calf    
%% Set variables for SortingHat code
% HbO2_run = [];
% Hb_run = [];
% Water_run = [];
% Lipid_run = [];
% Point_run = [];
Mua_run = [];
Mua_old = [];

waveSelect = 1;
patientVariable = i;
SortingHat_C
SortingHat_matlab
for str = 1:4
% for str = 1:7
    testLabels_C{str} = strcat(num2str(testwvs(1,str)),' C');
    testLabels_Matlab{str} = strcat(num2str(testwvs(1,str)),' M');
    % num2cell([testwvs(1,1:4) 'C'])
end

%% Analyze Metrics (Mua Diff, % Diff, Chro Diff, % Diff)
% Mua_run [C muaNew at testwvs (7) Matlab muaNew at testwvs (7)
% Mua_old [muaOld at testwvs(7) repeat]
clear Mua_Diff Mua_Perdiff Mua_STD;
Mua_Diff = abs(Mua_run - Mua_old);
% hold on
% plot(Mua_run(1,1:7),'r*')
% plot(Mua_run(1,8:14),'b*')

% plot(Mua_old)
Mua_Perdiff = 100*(Mua_Diff./Mua_old);
Mua_STD = std(Mua_Diff);
% figure;
figure('units','normalized','outerposition',[0 0 1 1])
% subplot(2,1,1);
% title('C Code')
% plot(repmat(testwvs(1,:),[length(Mua_Diff) 1]),Mua_Diff(:,1:7))%,'*r')
hold on
% boxplot(Mua_Diff(:,1:4),testwvs(1,1:4))

boxplot(Mua_Diff(:,1:8),[testLabels_C testLabels_Matlab]);%testwvs(1,1:4))
% boxplot(Mua_Diff(:,1:14),[testLabels_C testLabels_Matlab]);%testwvs(1,1:4))
title(strcat('Matlab vs C Code Residual {\mu}_a P',num2str(patientVariable)))
% boxplot(Mua_Diff(:,1:7),testwvs(1,1:7))
% title('C Code')
% boxplot(Mua_Diff(:,8:14),testwvs(1,1:7)+5)
% xlabel('Wavelength'),ylabel('{\mu}_a {\Delta}')
% subplot(2,1,2);
% plot(repmat(testwvs(1,:),[length(Mua_Diff) 1]),Mua_Diff(:,8:14),'^b')
% boxplot(Mua_Diff(:,8:14),testwvs(1,1:7))
% boxplot(Mua_Diff(:,5:8),testwvs(1,1:4))
% boxplot(Mua_Diff(:,1:14),[testwvs(1,1:7),testwvs(1,1:7)])
% title('Matlab Code')
% scatter(testwvs(1,:),Mua_Diff(1:7,:),50)
% scatter(testwvs(1,:),Mua_Diff(8:end,:),50)
xlabel('Wavelength'),ylabel('{\mu}_a {\Delta}') 
cd('C:\Users\Kyle\Documents\MATLAB\Analysis Pictures')
h = gcf;
savestring = strcat('Mua Difference ',num2str(patientVariable),tissue,'.tif');
saveas(h,savestring);
% close all
% figure;
figure('units','normalized','outerposition',[0 0 1 1])
% subplot(2,1,1);
hold on
% boxplot(Mua_Perdiff(:,1:7),testwvs(1,1:7))
% boxplot(Mua_Perdiff(:,1:4),testwvs(1,1:4))
% title('C Code')
% xlabel('Wavelength'),ylabel('{\mu}_a {\Delta}')
% subplot(2,1,2);
% plot(repmat(testwvs(1,:),[length(Mua_Diff) 1]),Mua_Diff(:,8:14),'^b')
% boxplot(Mua_Perdiff(:,8:14),testwvs(1,1:7))
% boxplot(Mua_Perdiff(:,5:8),testwvs(1,1:4))


boxplot(Mua_Perdiff(:,1:8),[testLabels_C testLabels_Matlab]);%testwvs(1,1:4))
% boxplot(Mua_Diff(:,1:14),[testLabels_C testLabels_Matlab]);%testwvs(1,1:4))


% boxplot(Mua_Diff(:,1:14),[testwvs(1,1:7),testwvs(1,1:7)])
% title('Matlab Code')
% title('Matlab vs C Code % Diff Residual {\mu}_a')
title(strcat('Matlab vs C Code % Diff Residual {\mu}_a P',num2str(patientVariable)))
% scatter(testwvs(1,:),Mua_Diff(1:7,:),50)
% scatter(testwvs(1,:),Mua_Diff(8:end,:),50)
xlabel('Wavelength'),ylabel('{\mu}_a {\Delta}') 
cd('C:\Users\Kyle\Documents\MATLAB\Analysis Pictures')
h = gcf;
savestring = strcat('Mua Perc Difference ',num2str(patientVariable),tissue,'.tif');
saveas(h,savestring);
% close all
% for j = 1:14
for j = 1:8
% fprintf(fid,'%5d \n',mean(Mua_Diff(j)));
fprintf(fid,'%15.12f \t',mean(Mua_Diff(:,j)));
fprintf(fid2,'%12.8f \t',mean(Mua_Perdiff(:,j)));
fprintf(fid3,'%12.8f \t',Mua_STD(j));
% fprintf(fid2,'%5d \n',mean(Mua_Perdiff(j)));
end
fprintf(fid,'\r\n');
fprintf(fid2,'\r\n');
fprintf(fid3,'\r\n');
end
fclose('all');

%%
% Abs values
HbO2_run =-abs(HbO2_run);
Hb_run =-abs(Hb_run);
Lipid_run =-abs(Lipid_run);
Water_run =-abs(Water_run);
%%
% Chromo Diff runs [mean min max std]
CHbO2 = HbO2_run(:,1:2:end);
CHb = Hb_run(:,1:2:end);
CLipid = Lipid_run(:,1:2:end);
CWater = Water_run(:,1:2:end);

HbO2 = HbO2_run(:,2:2:end);
Hb = Hb_run(:,2:2:end);
Lipid = Lipid_run(:,2:2:end);
Water = Water_run(:,2:2:end);
for p = 1
    
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
hold on
plot(-HbO2(p,:),'*r')
plot(-CHbO2(p,:),'^b')
title('HbO2')
hline = refline([0 0]);
set(hline,'Color','k')
xlabel('Patient'),ylabel('{\Delta} broadband')

subplot(2,2,2)
hold on
plot(-Hb(p,:),'*r')
plot(-CHb(p,:),'^b')
title('Hb')
hline = refline([0 0]);
set(hline,'Color','k')
xlabel('Patient'),ylabel('{\Delta} broadband')

subplot(2,2,3)
hold on
plot(-Lipid(p,:),'*r')
plot(-CLipid(p,:),'^b')
title('Lipid')
hline = refline([0 0]);
set(hline,'Color','k')
xlabel('Patient'),ylabel('{\Delta} broadband')

subplot(2,2,4)
hold on
plot(-Water(p,:),'*r')
plot(-CWater(p,:),'^b')
title('Water')
legend('Matlab', 'C Code')
hline = refline([0 0]);
set(hline,'Color','k')
xlabel('Patient'),ylabel('{\Delta} broadband')
end
