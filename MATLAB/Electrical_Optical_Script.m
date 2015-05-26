% Function: 
% Fusion of data from electrical and optical camera tracking
% Diffuse Optical Spectroscopic Imaging lab
% Kyle Cutler
% Beta v1 created 11/18/13
clear
clc
%% Read in data from Mouse

testname = 'MouseIMU_legsquare_1031.txt';
mouse_Data = load(testname);
mouse.rotational = mouse_Data(:,2:5);
mouse.displacement = mouse_Data(:,6:8);
mouse.time = mouse_Data(:,1);
clear mouse_Data testname
%% Read in data from Camera
clc
clear
meshname = 'exp_mesh4.ply';
[optical.vertex,optical.faces] = read_ply(meshname);
[d,c]=plyread2('exp_mesh4.ply')
%trisurf(optical.faces,optical.vertex(:,1),optical.vertex(:,2),optical.vertex(:,3));
%  %% Flatten 3D surface
faces = optical.faces;
vertex = optical.vertex;
name = 'lola';
options.name = name;
A = triangulation2adjacency(faces);
%%
options.boundary = 'circle';
options.laplacian = 'conformal';
options.method = 'parameterization';
options.verb = 0;
xy = compute_parameterization(vertex,faces,options);
%%
clf;
subplot(1,2,1);
plot_mesh(vertex,faces,options); shading faceted; axis tight;
subplot(1,2,2);
plot_graph(A,xy1,'k.-'); axis tight;
%%
options.method = 'freeboundary';
xy1 = compute_parameterization(vertex,faces,options);
%%
% display the parameterization
% clf;
% subplot(1,2,1);
% plot_graph(A,xy,'k.-'); axis tight;
% title('Fixed boundary');
% subplot(1,2,2);
plot_graph(A,xy1,'k.-'); axis tight;
title('Free boundary');


%% Select area/nodes for fiducial

%% Scale and align mouse path data

%% Project mouse path data on to flat surface based on fiducial

%% Unflatten 3D surface

%% Generate time,location,node data

%% Generate node normal vectors

%% Convert mouse quaternion to rotational matrix

%% Input starting orientation into rotational matrix

%% Comparison of two data sources
% 
%   for x = 1:10
%       disp(x)
%   end
% 