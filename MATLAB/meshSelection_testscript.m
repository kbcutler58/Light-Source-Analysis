clear
clc
load 3dDataTom.mat
vertices = meshBrt.vertices;
faces = meshBrt.faces;
clear meshBrt

bounds = [-50 50 -240 -150 -80 75];
[output1,output2] = meshSelection(vertices,faces,bounds);
% scatter3(output1(:,1),output1(:,2),output1(:,3))