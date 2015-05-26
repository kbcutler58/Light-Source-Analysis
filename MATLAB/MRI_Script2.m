%%Breast MRI Processing Script 9/6/13
% Kyle Cutler
%% Load data
clear
clc
load 3dDataTom.mat
vertices = meshBrt.vertices;
faces = meshBrt.faces;
clear meshBrt

% %Plot Original Data
trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3));
% trimesh(faces,vertices(:,1),vertices(:,2),vertices(:,3));
title('Original 3D Data'),xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')

%% Subselection of Data
% Bounds [lowerX upperX lowerY upperY lowerZ upperZ]
% Note that data in Y and Z switched in some sets
bounds = [-20 20 -220 -180 -20 20];
[subVertices,subFaces,interVertices,interFaces,inter2Faces,vertNumbers] = meshSelection(vertices,faces,bounds);

% %Plot Subselection
% trisurf(subFaces,subVertices(:,1),subVertices(:,2),subVertices(:,3));
% title('Subsampled Section'),xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')

%% Place measurement points with mouse
% Make sure to export cursor points at cursor_info
trisurf(subFaces,subVertices(:,1),subVertices(:,2),subVertices(:,3));
title('Subsampled Section'),xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')
%% Exp1

    for i=1:length(cursor_info)
        locationMatch=cursor_info(1,i).Position;
        a = (subVertices(:,1)==locationMatch(1));
        if sum(a) > 1
            b = subVertices(:,2)==locationMatch(2);
            if sum(b) > 1
                c = subVertices(:,3)==locationMatch(3);
                if sum(c) > 1
                    index(i)=max(find(c));
                else
                    index(i)=find(c);
                end
            else
                index(i)=find(b);
            end
        else
            index(i)=find(a);
        end
    end
    
    
    % To be taken out later when real data is available
    % Assign random values to measurement points
    values=rand(length(cursor_info),1);
    indexValues=zeros(length(subVertices),1);
    indexValues(index)=values;
    
    % Use measured points to make a nearest neighbor interpolation on mesh
    F = scatteredInterpolant(subVertices(index,1),subVertices(index,2),subVertices(index,3),indexValues(index),'nearest');
    
    % Apply interpolation to entire subsection of mesh
    Vq=F(subVertices(:,1),subVertices(:,2),subVertices(:,3));

%% Measurement and Scattered Interpolant Method
% F = scatteredInterpolant(x,y,z,v) creates 3D interpolant v = F(x,y,z)
% if exist('cursor_info','var') == 0;
%     load saved_cursor.mat
% else
% end
    % Match up mouse input with mesh elements
    for i=1:length(cursor_info)
        locationMatch=cursor_info(1,i).Position;
        a = (subVertices(:,1)==locationMatch(1));
        if sum(a) > 1
            b = subVertices(:,2)==locationMatch(2);
            if sum(b) > 1
                c = subVertices(:,3)==locationMatch(3);
                if sum(c) > 1
                    index(i)=max(find(c));
                else
                    index(i)=find(c);
                end
            else
                index(i)=find(b);
            end
        else
            index(i)=find(a);
        end
    end
    
    
    % To be taken out later when real data is available
    % Assign random values to measurement points
    values=10*rand(length(cursor_info),1);
    indexValues=zeros(length(subVertices),1);
    indexValues(index)=values;
    
    % Use measured points to make a nearest neighbor interpolation on mesh
    F = scatteredInterpolant(subVertices(index,1),subVertices(index,2),subVertices(index,3),indexValues(index),'nearest');
    
    % Apply interpolation to entire subsection of mesh
    Vq=F(subVertices(:,1),subVertices(:,2),subVertices(:,3));
% %Plot measurement points
% close all
% trisurf(subFaces,subVertices(:,1),subVertices(:,2),subVertices(:,3),indexValues);
% title('Subsampled Section with point measurements'),xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')

% %Plot interpolant of subsection
% figure
% trisurf(subFaces,subVertices(:,1),subVertices(:,2),subVertices(:,3),Vq);
% title('Subsampled Section with expanded area measurements'),xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')

%Plot subsection and original dataset
clear C
C=zeros(length(vertices),1);
for i=1:length(vertices)
    if vertices(i,2) <= -140
        if vertices(i,3) >= -10 && vertices(i,3) <= 15
            if vertices(i,1) >= 10 && vertices(i,1) <= 25
                C(i,1)= 0;
            else
            end
        else
        end
    else
    C(i,1)= 0;
    end    
end
%%
figure;
hold on
trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3),C)
trisurf(subFaces,subVertices(:,1),subVertices(:,2),subVertices(:,3),Vq);
title('Subsampled Section overlay onto full dataset'),xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')
hold off

figure;
hold on
trimesh(faces,vertices(:,1),vertices(:,2),vertices(:,3),C)
trimesh(subFaces,subVertices(:,1),subVertices(:,2),subVertices(:,3),Vq);
title('Subsampled Section overlay onto full dataset'),xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')
hold off

figure;
hold on 
trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3),C)
trimesh(subFaces,subVertices(:,1),subVertices(:,2),subVertices(:,3),Vq);
title('Subsampled Section overlay onto full dataset'),xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')
hold off
%%
figure;
hold on
trimesh(faces,vertices(:,1),vertices(:,2),vertices(:,3),C)
trisurf(subFaces,subVertices(:,1),subVertices(:,2),subVertices(:,3),Vq);
title('Subsampled Section overlay onto full dataset'),xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')
hold off
%% Experimental

%{ 

%Everything greater than -140 in Y = 0
%TOI ranges from 0 - 10 with highest values at tumor

% clear xver yver zver
% subsampleFactor=5;
% j=1;
% for i=1:subsampleFactor:length(vertices)
%     xver(j,1)=vertices(i,1);
%     yver(j,1)=vertices(i,2);
%     zver(j,1)=vertices(i,3); 
%     j=j+1;
% end
% for i=1:subsampleFactor:length(faces)
%     newfaces(j,:)=faces(i,:);
%     j=j+1;
% end
% 
% 
% Vq = F(xver,yver,zver);
% % mesh(xver,yver,zver,Vq);
% triver = delaunay(xver,yver,zver);
% trisurf(triver,xver,yver,zver,Vq);

% %F = scatterInterpolant(vertices(:,1),vertices(:,2),vertices(:,3),vertices(:,3)*10);
% %VI = interp3(values,vertices(:,1),vertices(:,2),vertices(:,3));
% C=vertices(:,2)/min(vertices(:,2));
% F = TriScatteredInterp(vertices(:,1),vertices(:,2),vertices(:,3),C);
% 
% trimesh(faces,vertices(:,1),vertices(:,2),vertices(:,3));
% figure;
% trimesh(faces,vertices(:,1),vertices(:,2),vertices(:,3),F.V);
%}