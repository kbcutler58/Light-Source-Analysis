% plot(time2,perfusion2,'y');
% hold on
% [peaks,location]=findpeaks(perfusion2,'minpeakdistance',15);
% scatter(time2(location),perfusion2(location));
% [peaks2,location2]=findpeaks(perfusion2.^-1,'minpeakdistance',15);
% scatter(time2(location2),perfusion2(location2));
%% Selection of time points
close all
timepoint = 80; %time threshold
[newtime_i]=find(time2>timepoint);
time_set = time2(newtime_i);
perf_set = perfusion2(newtime_i);
% plot(time_set,perf_set)
%% Selction of laserdata
close all
hold on
laserpoint = 100; %perfusion threshold
[newtime_i2]=find(perf_set>laserpoint);
time_set = time_set(newtime_i2);
perf_set = perf_set(newtime_i2);
plot(time_set,perf_set)
[peaks,location]=findpeaks(perf_set,'minpeakdistance',15);
[peaks2,location2]=findpeaks(perf_set.^-1,'minpeakdistance',15);
scatter(time_set(location),perf_set(location),'y')
scatter(time_set(location2),perf_set(location2),'k')
%%

%location = peak index
%location2= trough index
newpoints = zeros(length(location)-1,1);
newpoints2 = zeros(length(location)-1,1);
for i=1:length(location)-1
    u = location(i);
    v = location(i+1);
    [z,location3]=find(location2 < v & location2 > u, 1, 'last');
%     if isempty(z)   
    newpoints(i)= peaks(i)- perf_set(location2(z));
    newpoints2(i) = (peaks(i)+perf_set(location2(z)))/2;
end
hold on
plot(time_set(location(1:end-1)),newpoints)
plot(time_set(location(1:end-1)),newpoints2)
% plot(newpoints2)

    