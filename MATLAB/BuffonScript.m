% MC Buffon Simulation
% BME 230B HW3 A
% Kyle Cutler

%%Generate Random Variables between 0 and 1
clear a
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));
a=rand(100000,3);
%%Create Sticks
x = a(:,1);
y = a(:,2);
direction=a(:,3)*2*pi;
long =.35;
[ends] = stickgen(x,y,direction,long);

% plot(x,y,'r*',ends(:,1:2),ends(:,3:4),'*b'), xlim([0;1]), ylim([0;1]) 
% for j=1:length(ends)
%     hold on
%     plot(ends(j,1:2),ends(j,3:4)), xlim([0;1]), ylim([0;1])
% end
% 
%  for j=1:length(ends)
%     hold on
%     plot(ends(j,1:2)-.5,ends(j,3:4)), xlim([-.5;.5]), ylim([0;1])
% end
% plot([0 0],[0 1],'r')

p = stickcheck(ends);
pie = (2*long)/(p*1)