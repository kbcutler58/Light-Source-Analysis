x_raw = data(:,2); % raw x data
y_raw = data(:,3); % raw y data
yaw = deg2rad(data(:,4)); % convert yaw to rads
yaw2 = yaw - yaw(1); % offset angle to zero
x_location2 = 0; % initial x with no offset
y_location2 = 0; %initial y with no offset

for i=1:length(x_raw)
    % Calculate rotation matrix 
    correction_factor = [ cos(yaw2(i)) -sin(yaw2(i)); sin(yaw2(i)) cos(yaw2(i))];
    % Matrix for raw x and y
    path_matrix = [x_raw(i),y_raw(i)];
    % Multiply x,y matrix and rotation matrix
    path_temp = path_matrix*correction_factor;
    % Add temp path to previous location
    corrected_path2(i,:) = [x_location2 + path_temp(1) , y_location2 + path_temp(2)];
    % Update Location of x and y 
    x_location2 = corrected_path2(i,1);
    y_location2 = corrected_path2(i,2);
    scatter(corrected_path2(1:i,1),corrected_path2(1:i,2),100,-colordata(1:i,1))
    axis equal
    colorbar
    title('2D Path Visualization')
    F(i) = getframe(gcf);
end
%%
for i=1:length(x_raw)
    scatter(corrected_path2(1:i,1),corrected_path2(1:i,2),150,-colordata(1:i,1))
    axis equal
    colorbar
    title('2D Path CW Visualization')
    F(i) = getframe(gcf);
end    