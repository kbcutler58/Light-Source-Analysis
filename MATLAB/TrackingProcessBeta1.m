clear
clc
data = load('mousedata1.mat');
data_10mm = data.data10mm2;
data_5mm = data.data5mm2;
clear data
%%
x1 = data_5mm(:,2);
for i=1:length(x1)
    if x1(i)>129
        x1(i) = x1(i)-256;
    else
    end
end
x2 = data_10mm(:,2);
for i=1:length(x2)
    if x2(i)>129
        x2(i) = x2(i)-256;
    else
    end
end
%%
test1 = trapz(x1);
test2 = trapz(x2);

%%
x3 = skim10cm1(:,2);
x4 = skim5cm1(:,2);
for i=1:length(x3)
    if x3(i)>129
        x3(i) = x3(i)-256;
    else
    end
end
for i=1:length(x4)
    if x4(i)>129
        x4(i) = x4(i)-256;
    else
    end
end
test3 = trapz(x3);
test4 = trapz(x4);