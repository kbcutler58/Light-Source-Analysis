% % Calculate Chromophore Concentration Matrix
clear TotalChrom
% for wvcase = 0:28
%     waverange = (1:29:length(HbO2_run))+wvcase;
for wvcase = 0:16
    waverange = (1:17:length(HbO2_run))+wvcase;
TotalChrom(1,wvcase+1) = sum((HbO2_run(1,waverange).*Point_run),2)./sum(Point_run);
TotalChrom(2,wvcase+1) = sum((HbO2_run(4,waverange).*Point_run),2)./sum(Point_run);
TotalChrom(3,wvcase+1) = min(HbO2_run(2,waverange));
TotalChrom(4,wvcase+1) = max(HbO2_run(3,waverange));

TotalChrom(5,wvcase+1) = sum((Hb_run(1,waverange).*Point_run),2)./sum(Point_run);
TotalChrom(6,wvcase+1) = sum((Hb_run(4,waverange).*Point_run),2)./sum(Point_run);
TotalChrom(7,wvcase+1) = min(Hb_run(2,waverange));
TotalChrom(8,wvcase+1) = max(Hb_run(3,waverange));

TotalChrom(9,wvcase+1) = sum((Water_run(1,waverange).*Point_run),2)./sum(Point_run);
TotalChrom(10,wvcase+1) = sum((Water_run(4,waverange).*Point_run),2)./sum(Point_run);
TotalChrom(11,wvcase+1) = min(Water_run(2,waverange));
TotalChrom(12,wvcase+1) = max(Water_run(3,waverange));

TotalChrom(13,wvcase+1) = sum((Lipid_run(1,waverange).*Point_run),2)./sum(Point_run);
TotalChrom(14,wvcase+1) = sum((Lipid_run(4,waverange).*Point_run),2)./sum(Point_run);
TotalChrom(15,wvcase+1) = min(Lipid_run(2,waverange));
TotalChrom(16,wvcase+1) = max(Lipid_run(3,waverange));
end