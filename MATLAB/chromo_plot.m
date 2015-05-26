

totalPoints = size(fieldnames(selectedDataFinal.(fields3{1}).(fields4{1})),1);
% totalPoints = length(chromDif)/29
disp(totalPoints);

clear Hb HbO2 Water Lipid
 
titlestring = strcat('Patient ',num2str(patientSelect(1)),' Case ',num2str(waveSelect(1)),':',num2str(waveSelect(end)),' Conc Changes');

% figure('units','normalized','outerposition',[0 0 1 1]);
% pointrange(:,1) = totalPoints*(waveSelect-1)+1;
% pointrange(:,2) = totalPoints+totalPoints*(waveSelect-1);
% wvcombo = [];
% for i = 1:length(waveSelect)
%     wvcombo = [wvcombo; repmat(i,pointrange(i,2)-pointrange(i,1)+1,1)];
% end
% 
% boxplot(chromDif(:,3),wvcombo)

% scatter(wvcombo,chromDif(:,1))
subplot(2,2,1)
for i = 1:length(waveSelect)
    hold on
    wvcombo = repmat(i,totalPoints,1 );
    scatter(wvcombo,chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),1))
%     scatter(wvcombo,chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),1))
    title(strcat(titlestring,' HbO2'))
    xlabel('Wavelength Test Case')
    ylabel('HbO2 difference in {\mu}M')
    ax = gca;
    set(ax,'YLim',[(-max(abs(chromDif(:,1))))-1 max(abs(chromDif(:,1)))+1])
    disp(strcat('Case ',num2str(i,'%02d')));
    HbO2(1,i) = mean(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),1));
    HbO2(2,i) = min(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),1));
    HbO2(3,i) = max(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),1));
    HbO2(4,i) = std(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),1));
%     HbO2(5,i) = totalPoints;
end
    hline = refline([0 0]);
    set(hline,'Color','k')
    
subplot(2,2,2);
for i = 1:length(waveSelect)
    hold on
    wvcombo = repmat(i,totalPoints,1 );
    scatter(wvcombo,chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),2))
    title(strcat(titlestring,' Hb'))
    xlabel('Wavelength Test Case')
    ylabel('Hb difference in {\mu}M')
    ax = gca;
    set(ax,'YLim',[(-max(abs(chromDif(:,2))))-1 max(abs(chromDif(:,2)))+1])
    Hb(1,i) = mean(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),2));
    Hb(2,i) = min(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),2));
    Hb(3,i) = max(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),2));
    Hb(4,i) = std(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),2));
%     Hb(5,i) = totalPoints;
end
    hline = refline([0 0]);
    set(hline,'Color','k')
    
subplot(2,2,3)
for i = 1:length(waveSelect)
    hold on
    wvcombo = repmat(i,totalPoints,1 );
    scatter(wvcombo,chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),3))
    title(strcat(titlestring,' Water Frac'))
    xlabel('Wavelength Test Case')
    ylabel('Water Frac difference')
    ax = gca;
    set(ax,'YLim',[(-max(abs(chromDif(:,3))))-1 max(abs(chromDif(:,3)))+1])
    Water(1,i) = mean(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),3));
    Water(2,i) = min(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),3));
    Water(3,i) = max(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),3));
    Water(4,i) = std(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),3));
%     Water(5,i) = totalPoints;
end

    hline = refline([0 0]);
    set(hline,'Color','k')
    
subplot(2,2,4)
for i = 1:length(waveSelect)
    hold on
    wvcombo = repmat(i,totalPoints,1 );
    scatter(wvcombo,chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),4))
    title(strcat(titlestring,' Lipid Frac'))
    xlabel('Wavelength Test Case')
    ylabel('Lipid Frac difference')
    ax = gca;
    set(ax,'YLim',[(-max(abs(chromDif(:,4))))-1 max(abs(chromDif(:,4)))+1])
    Lipid(1,i) = mean(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),4));
    Lipid(2,i) = min(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),4));
    Lipid(3,i) = max(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),4));
    Lipid(4,i) = std(chromDif(totalPoints*(i-1)+1:totalPoints+totalPoints*(i-1),4));
%     Lipid(5,i) = totalPoints;
end
    hline = refline([0 0]);
    set(hline,'Color','k')
 
% cd('C:\Users\Kyle\Downloads\DOSI data for wv LG testing\Analysis Pictures')
% h = gcf;
% savestring = strcat('ChromDist ',num2str(patientSelect),'case',num2str(waveSelect(1)),' Calf','.jpg');
% savestring = strcat('ChromDist ',num2str(patientSelect),' Case All','Tumor Breast 4wv','.jpg');
% saveas(h,savestring);
% savestring = strcat('ChromDist ',num2str(patientSelect),' Case All','Tumor Breast 4wv','.fig');
% saveas(h,savestring);

