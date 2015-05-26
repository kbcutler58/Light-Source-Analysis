

%% Subselecting VC

% Only use z range VC
% VC1 = vertex_correspondence(:,vertselect);
VC1 = vertex_correspondence;

% Scale VC in x and y
ScaleFactor = 4.5;
VC2 = VC1 * ScaleFactor;
U = -60;
V = -5;
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

%% Sort and index VC
VC_final_only = VC_final(:,VC_mask);
VC_final_x = VC_final(1,VC_mask);
VC_final_y = VC_final(2,VC_mask);

[VC_sorted_x,index_x] = sort(VC_final_x,2);
[VC_sorted_y,index_y] = sort(VC_final_y,2);

%% Tag sorted VC with 2D locations
y_range3 = fliplr(y_range);
%%
dx = 2.5/4;
dy = 2.5/4;
x_range2= min(x_range):dx:max(x_range);
y_range2= min(y_range3):dy:max(y_range3);
%%
% x_index = zeros(length(x_range),1);
% y_index = zeros(length(y_range),1);
% x_index = 1;
% y_index = 1;
% x_index(length(x_range))=length(index_x);
% y_index(length(y_range))=length(index_y);
% 
% for i=1:length(x_range)-1
%     temp = find(VC_sorted_x>(x_range(i+1)),1);
%     if isempty(temp)
%         temp = length(index_x);
%     end
%     x_index(i+1) = temp;
%     index_x(2,x_index(i):x_index(i+1))=i;
%     
% end
% 
% for j=1:length(y_range)-1
%     temp = find(VC_sorted_y>(y_range(j+1)),1);
%     if isempty(temp)
%         temp = length(index_y);
%     end
%     y_index(j+1) = temp;
%     index_y(2,y_index(j):y_index(j+1))=j;
% end

%%
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

%% Color pixels according to 
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

%%
trimesh(faces,vert(1,:),vert(2,:),vert(3,:),colorData2);
set(gca,'clim',colorbaroriginal)
colorbar
view(2)