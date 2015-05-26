% Create animation of CW/Tracking plot

for i=1:length(data_x)
    scatter(corrected_path_offset(1:i,1),corrected_path_offset(1:i,2),100,-peak_power(1:i,3))
    axis equal
    colorbar
    title('2D Path Visualization 20khz peak')
    Movie_11khz(i) = getframe(gcf);
    close all
    
    scatter(corrected_path_offset(1:i,1),corrected_path_offset(1:i,2),100,-peak_power(1:i,1))
    axis equal
    colorbar
    title('2D Path Visualization 11khz peak')
    Movie_20khz(i) = getframe(gcf);
    close all
end

%%
% movie2avi(F,'CW_Tracking2.avi','fps','150')