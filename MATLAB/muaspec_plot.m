% Plot of Discrete Mua fit vs Absorption and Broadband fit


% measure_type = ' Breast';
% measure_type = ' Breast Tumor';
% measure_type = ' Abdomen';
measure_type = ' Calf';
% measure_type = ' Bicep'

% Options to plot just one
range = 1;

% range = 1:45;

% % Options to automatically create full screen plots
% figure('units','normalized','outerposition',[0 0 1 1]);

% Option to plot all mua spectrum
% range = 1:size(muaOld,1);

% domain = find(muaOld(1,:) ~= 0);
% titlestring = strcat('Patient ',num2str(patientSelect(1)),':',num2str(patientSelect(end)),' Wavelength Case ',num2str(waveSelect(1)),':',num2str(waveSelect(end)),' Broadband vs Discrete Point Fit');
titlestring = strcat('Patient ',num2str(patientSelect(1)),' Wavelength Case ',num2str(waveSelect(1)),' Broadband vs Discrete Point Fit ',measure_type);

% for r= 1:length(range);
% figure;
r = range;
hold on
domain = find(muaOld(r,:) ~= 0);
plot(wvs(domain),muaNew(r,domain),'b',wvs(domain),muaOld(r,domain),'r',wvs(domain),muaOldFit(r,domain),'g','LineWidth',1.5),
[a,b,c] = intersect(testwvs(1,:),wvs);
scatter(testwvs(1,:),muaOld(r,c),50,[1 0 0],'fill')
hold off
% end
title(titlestring)
xlabel('wavelength (nm)'), ylabel('{\mu}_a (1/mm)');
leg1 = legend('Discrete Chromophore Fit','Absorption Spectrum','Broadband Chromophore Fit','Discrete Absorption Values');
ax = gca;
set(leg1,'Location','Best')


% cd('C:\Users\Kyle\Downloads\DOSI data for wv LG testing\Analysis Pictures')
% h = gcf;
% % savestring = strcat('ChromDist ',num2str(patientSelect),'case',num2str(waveSelect(1)),' Calf','.jpg');
% savestring = strcat('Sample Spectrum Patient',num2str(patientSelect),' Point ',num2str(range),' Case 1','Tumor Breast 4wv','.jpg');
% saveas(h,savestring);
% savestring = strcat('Sample Spectrum Patient',num2str(patientSelect),' Point ',num2str(range),' Case 1','Tumor Breast 4wv','.fig');
% saveas(h,savestring);

% close all