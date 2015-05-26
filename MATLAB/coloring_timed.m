colorData2 = zeros(length(vert),1);
for i=1:length(x_range2)-1
    for j=1:length(y_range2)-1
        xset = index_x(1,x_index(i):x_index(i+1));
        yset = index_y(1,y_index(j):y_index(j+1));
        combinationSet = intersect(xset,yset);
        combocombo = VC_final_only(3,combinationSet);
        colorData2(combocombo) = colorData(241-j,i);
    end
end

trimesh(faces,vert(1,:),vert(2,:),vert(3,:),colorData2);
set(gca,'clim',colorbaroriginal)
colorbar
view(2)