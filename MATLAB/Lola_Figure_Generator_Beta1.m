%{
    DOSI Lab Lola 3D Figure Generator Script
        Requires stlread.m
            [faces,vertices,normals]=stlread(stlfile.stl)
        Requires cameraMatrix.m
            [transformation matrix] = cameraMatrix(4xM real coordinates,2xM
            2D coordinates)
    Inputs:
        Path to folder with figure files
        STL file of lola
        Vertex correspondence data file
    Outputs:
        3D figures    

    File Structure
        Load in STL file
        Load in data figure(s)
        Load/Create vertex correspondence
        3D to 2D transform
        Subselect vertex map
        Tag sorted vertex map with 2D locations
        Color vertices according to tag
        Display 3D Figure
        Save File
%}
clear

% Visualization Options
ScaleFactor = 4.5; %Size scaling
U = -60; % Left and Right translation
V = -5;  % Up and Down translation

%% Load in STL file

stl_file = 'stl_test_file.stl';
[faces,vert,normals]=stlread(stl_file);
x = 100*vert(:,1);
y = 100*-vert(:,2);
z = 100-vert(:,3);
vert = [x,y,z];
% trimesh(faces,x,y,z);
% view(2)

%% Load in data figure

folder = 'C:\Users\Kyle\Desktop\MRI\8714-01-DOSI\Baseline_DOSI Images_scaled\figure files';
cd(folder);
folderlabel = dir('*.fig');
[pathstr,name,ext]=fileparts(folderlabel(1).name);
imagelabel = name;
pathtoimage = strcat(imagelabel, ext);
open(pathtoimage);
h = gcf;
[x_range,y_range,colorData]=getimage(h);
colorbaroriginal=get(gca,'clim');
close all


%% Create/Load correspondence
folder = 'C:\Users\Kyle\Downloads';
cd(folder);

datatip = load('cursor_info_oct6.mat');
cursor_info = datatip.cursor_info;
clear datatip

% %% Create vertex correspondence
% trimesh(faces,x,y,z);
% view(2);

% Select Points on Model
for i=1:length(cursor_info)
    transform1(1:3,i) =cursor_info(1,i).Position;
end
transform1(4,:) = 1;

% Assign 2D Correspondence

% 2xM matrix [ x1 x2 x3...xn,y1 y2 y3...yn]
% Original 2D Points (0,0) (0,10) (10,10) (10,0)
transform2=[0 0 10 10; 0 10 10 0];

% Generate Transformation Matrix

[correspondence_matrix] = cameraMatrix(transform1,transform2);

% Subselect area of mesh to map data (R or L)

zrange = 99.32;
vertselect=find(vert(:,3)>zrange);
vertcopy=zeros(length(vert),3);
vertcopy = vertcopy + 1000;
vertcopy(:,3)=1;
vertcopy(vertselect,:)= vert(vertselect,:);
vertcopy=vertcopy';
vertcopy(4,vertselect)=1;
vert=vert';
vert(4,:)=1;

%% 3D to 2D transform

% %Use this line instead of next if you do not want Z cutoff
% new_correspondence = correspondence_matrix*vert;

new_correspondence = correspondence_matrix*vertcopy;
vertex_correspondence = zeros(3,length(new_correspondence));

for i = 1:length(new_correspondence);
    if (new_correspondence(3,i) ~= 0)
        vertex_correspondence (1,i) = new_correspondence(1,i)/new_correspondence(3,i); 
        vertex_correspondence (2,i) = new_correspondence(2,i)/new_correspondence(3,i); 
        vertex_correspondence (3,i) = new_correspondence(3,i)/new_correspondence(3,i); 
    else
    end
end

Figure_Gen_exp

%% Subselecting VC

% Only use z range VC
% VC1 = vertex_correspondence(:,vertselect);
VC1 = vertex_correspondence;

% Scale VC in x and y
VC2 = VC1 * ScaleFactor;
VC2(1,:) = VC2(1,:) + U;
VC2(2,:) = VC2(2,:) + V;

% Select vertex correspondence that fits image
lowVC_xlimit = min(x_range);
highVC_xlimit = max(x_range);
lowVC_ylimit = min(y_range);
highVC_ylimit = max(y_range);

lowerVCx = (VC2(1,:)<lowVC_xlimit);
upperVCx = (VC2(1,:)>highVC_xlimit);
lowerVCy = (VC2(2,:)<lowVC_ylimit);
upperVCy = (VC2(2,:)>highVC_ylimit);

VC_final = zeros(2,length(VC2));
VC_mask1 = ((lowerVCx | upperVCx));
VC_mask2 = ((lowerVCy | upperVCy));
VC_mask = ((~VC_mask1) & (~VC_mask2));
VC_final(:,VC_mask) = VC2(1:2,VC_mask);
VC_final(3,:) = 1:length(VC_final);

% Sort and index VC
VC_final_only = VC_final(:,VC_mask);
VC_final_x = VC_final(1,VC_mask);
VC_final_y = VC_final(2,VC_mask);

[VC_sorted_x,index_x] = sort(VC_final_x,2);
[VC_sorted_y,index_y] = sort(VC_final_y,2);

%% Tag sorted VC with 2D locations
y_range3 = fliplr(y_range);
dx = 2.5/4;
dy = 2.5/4;
x_range2= min(x_range):dx:max(x_range);
y_range2= min(y_range3):dy:max(y_range3);

x_index = zeros(length(x_range2),1);
y_index = zeros(length(y_range2),1);
x_index = 1;
y_index = 1;
x_index(length(x_range2))=length(index_x);
y_index(length(y_range2))=length(index_y);

for i=1:length(x_range2)-1
    temp = find(VC_sorted_x>(x_range2(i+1)),1);
    if isempty(temp)
        temp = length(index_x);
    end
    x_index(i+1) = temp;
    index_x(2,x_index(i):x_index(i+1))=i;
    
end

for j=1:length(y_range2)-1
    temp = find(VC_sorted_y>(y_range2(j+1)),1);
    if isempty(temp)
        temp = length(index_y);
    end
    y_index(j+1) = temp;
    index_y(2,y_index(j):y_index(j+1))=j;
end

%% Color vertices according to Tag
colorData2 = zeros(length(vert),1);
for i=1:length(x_range2)-1
    for j=1:length(y_range2)-1
        xset = index_x(1,x_index(i):x_index(i+1));
        yset = index_y(1,y_index(j):y_index(j+1));
        combinationSet = intersect(xset,yset);
        combocombo = VC_final_only(3,combinationSet);
        colorData2(combocombo) = colorData(241-j,289-i);
    end
end

%% Display Figure
trimesh(faces,vert(1,:),vert(2,:),vert(3,:),colorData2);
set(gca,'clim',colorbaroriginal)
colorbar
view(2)

% Previous Method
%{ 
% %% Color vertices according to data figure 
% % Method 1
% 
% U = 0;
% V = 0;
% % vertex_correspondence = vertex_correspondence';
% vertex_correspondence2=vertex_correspondence(:,vertselect);
% vc = vertex_correspondence2'*10; %Scaled to mm from cm
% dx = 2.5/4;  %pixels per x distance
% dy = -2.5/4; %pixels per y distance
% Color_Vector = zeros(length(vert),1); %Vertex color vector
% Color_List2 = zeros(length(vert),1); %temp color vector
% testcase4=zeros(length(vert),1);
% for i=1:length(x_range)-1
%     for j=1:length(y_range)-1
%         for k=1:4
%             for l=1:4
%                 x1=x_range(i)+dx*(k-1)+U;
%                 x2=x_range(i)+dx*(k)+U;
%                 y1=y_range(j)+dy*(l-1)+V;
%                 y2=y_range(j)+dy*(l)+V;
%                 testcase1 = (vc(:,1) > x1) & (vc(:,1) < x2);
%                 testcase2 = (vc(:,2) > y2) & (vc(:,2) < y1);
%                 testcase3 = testcase1 & testcase2; 
%                 Color_List = testcase3 * colorData((j-1)*4+l,(i-1)*4+k);
%                 Color_List2(vertselect)=Color_List;
%                 testcase4(vertselect)=testcase3;
%                 if sum(Color_List) > 0
%                 Color_Vector(find(testcase4==1)) = Color_List2(find(testcase4==1));
%                 else
%                 end
%             end
%         end
%     end
% end
% 
% % Method 2
% % Precalculate boundary numbers that correspond to pixels
% % Add vertex index to VC
% % Sort VC by VC(x)
% % Sort VC by VC(y)
% % for x
% % for y
% % index1(i) = find (VC(x) == boundary)
% % index2(i) = find (VC(y) == boundary)
% % pixel tag = VC(index1(i)-index(i-1) & index2(i)-index2(i-1))
% % end
% % end
% 
% %% Display 3D Figure
% trimesh(faces,vert(:,1),vert(:,2),vert(:,3),Color_Vector);
% set(gca,'clim',colorbaroriginal)
% colorbar
% view(2)
% %% Save File

%}


