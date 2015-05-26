%% Tissue Matrix Setup
close all
clear
clc
Tissue_Mat=zeros(12,22);
a=...
[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1;
 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0;
 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1;
 1, 1, 1, 1, 2, 0, 0, 2, 2, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 0, 0, 0, 2, 2, 2, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1;
 1, 1, 1, 1, 0, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1;
 1, 1, 1, 1, 0, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1;
 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1;
 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
Tissue_Mat(2:11,2:21)=a;
clear a
% imagesc(Tissue_Mat)
%% Concentration Matrix and Diffusion Matrix
%{
close all
D=zeros(12,22);
for j = 2:21
    for i = 2:11
        if Tissue_Mat(i,j) == 0
            if Tissue_Mat(i+1,j) == 0
                temp(1) = 10;
            else
                temp(1) = 2.5;
            end
            if Tissue_Mat(i-1,j) == 0
                temp(2) = 10;
            else
                temp(2) = 2.5;
            end
            if Tissue_Mat(i,j+1) == 0
                temp(3) = 10;
            else
                temp(3) = 2.5;
            end
            if Tissue_Mat(i,j-1) == 0
                temp(4) = 10;
            else
                temp(4) = 2.5;
            end
            D(i,j) = mean(temp);
        else
            D(i,j) = 2.5;
        end
    end
end
D=D/100;
% Dstar=D(2:11,2:21);
%}
%Average D calculation
%Initial Concentration Matrix
u=zeros(12,22);
%BC's for Concentration Matrix Top Left =2 Bottom Right = 0
u(1,:)=2e-6;
u(12,:)=0;
u(:,1)=2e-6;
u(:,22)=0;
%Diffusion 3D Matrix
D = zeros(12,22,4);
for i = 2:11
    for j = 2:21
        if Tissue_Mat(i,j) == 0 && Tissue_Mat(i-1,j) == 0
            D(i,j,1) = 10;
        elseif Tissue_Mat(i,j) == 0 && Tissue_Mat(i+1,j) == 0
            D(i,j,2) = 10;
        elseif Tissue_Mat(i,j) == 0 && Tissue_Mat(i,j-1) == 0
            D(i,j,3) = 10;
        elseif Tissue_Mat(i,j) == 0 && Tissue_Mat(i,j+1) == 0
            D(i,j,4) = 10;
        end
        for k = 1:4
            if D(i,j,k) == 0
                D(i,j,k) = 2.5;
            else
            end
        end
    end
end
D=D*1e-3;
%{ 
Test Code for Diffusion Matrix
Tissue_Mat(2,2)
d(3,2,1)
Tissue_Mat(4,2)
d(3,2,2)
Tissue_Mat(1,2)
d(3,2,3)
Tissue_Mat(2,2)
d(3,2,4)
%}        
 %Dilution Factor
%{
% FTCS Diffusion without Metabolism
t_length=10000;
dudt=zeros(12,22);
t = 1;
concdif = 1;
while concdif > 1e-9 % 1:t_length
% for t = 1:t_length
    prevconcentration=sum(sum(u));
    for i = 2:11%:11
        for j =2:21%:21
            dudx1 = D(i,j,1)* (u(i-1,j));
            dudx2 = D(i,j,2)* (u(i+1,j));
            dudy1 = D(i,j,3)* (u(i,j-1));
            dudy2 = D(i,j,4)* (u(i,j+1));
            dudt(i,j) = (dudx1+dudx2+dudy1+dudy2+(-u(i,j))*4*mean(D(i,j,:)));
        end
    end
    
    for i = 2:11
        for j = 2:21
            u(i,j) = u(i,j)+ .5 * dudt(i,j);%(.5 * D(i,j) * dudt(i,j));
        end
    end
    u(1,:)=2e-6;
    u(12,:)=0;
    u(:,1)=2e-6;
    u(:,22)=0;
    
    newconcentration=sum(sum(u));
    concdif(t) = newconcentration - prevconcentration;
    t=t+1;
end
%}
%% Concentration Differences over time (Convergence)

close all
figure;
subplot(2,1,1)
imagesc(u),title('Concentration of Oxygen'), xlabel('tissue width in microns'), ylabel('tissue height in microns')%concentration(:,:,x)), title('Concentration of Oxygen')
subplot(2,1,2)
imagesc(Tissue_Mat), title('Tissue Type')

%% FTCS Steady State Solver with Cell Metabolism

t_length = 15000;
% t = 1;
concdif = 1;
%while concdif > 1e-9 % 1:t_length
for t = 1:t_length
    prevconcentration=sum(sum(u));
    for i = 2:11%:11
        for j =2:21%:21
            dudx1 = D(i,j,1)* (u(i-1,j));
            dudx2 = D(i,j,2)* (u(i+1,j));
            dudy1 = D(i,j,3)* (u(i,j-1));
            dudy2 = D(i,j,4)* (u(i,j+1));
            dudt(i,j) = (dudx1+dudx2+dudy1+dudy2+(-u(i,j))*4*mean(D(i,j,:)));
        end
    end
    
    for i = 2:11
        for j = 2:21
            u(i,j) = u(i,j)+ .5 * dudt(i,j);%(.5 * D(i,j) * dudt(i,j));
            if u(i,j) > 2e-9  && Tissue_Mat(i,j) == 2
                u(i,j) = u(i,j) - 1e-9;
            end
        end
    end
    
    u(1,:)=2e-6;
    u(12,:)=0;
    u(:,1)=2e-6;
    u(:,22)=0;
    
    newconcentration=sum(sum(u));
    concdif(t) = newconcentration - prevconcentration;
    %t=t+1;
end

%% Convergence
close all
figure;
plot(concdif(2:length(concdif))), title('Convergence of solution'), xlabel('time steps'), ylabel('difference in sum of concentrations between time steps');
%%
figure;
subplot(2,1,1)
imagesc(u),title('Concentration of Oxygen')%concentration(:,:,x)), title('Concentration of Oxygen')
subplot(2,1,2)
imagesc(Tissue_Mat), title('Tissue Type')


