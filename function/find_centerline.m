function [x_axis,y_axis,max_velocity] = find_centerline(velocity_magnitude, zoom_scale)

% find the trajectory of the jets and output the location of tracjectory

velocity_result = imresize(velocity_magnitude, zoom_scale);      % resize the image
[row, ~] = size(velocity_result);
velocity_result(velocity_result<1e-3) = 0;       %low pass the small value

find_edge = edge(velocity_result,'Canny');       %find the edge and the range of jets
find_null = find(sum(find_edge) == 0);
find_edge(:,find_null) = [];                     %delete the data out of the jets
velocity_result(:,find_null) = [];
[~, column] = size(find_edge);
index_v = zeros(1,column);
max_velocity = zeros(1,column);
for j = 1 : column
    [max_velocity(:,j), index] = max(velocity_result(:,j));
    index_v(:,j) = index;                        %record the location of maximum velocity
    edge_index = find(find_edge(:,j) == 1);
    edge_index = [min(edge_index); max(edge_index)];       %the horizontal range of jets 
end

% zoom out the max velocity to the original data
max_velocity = imresize (max_velocity, 1 / zoom_scale);

% zoom out the index to the original data (size and value)
x_axis = (1 : column / zoom_scale)';
index_v_zoomout = imresize (index_v(:), 1 / zoom_scale);
y_axis = (row -  index_v_zoomout) ./ zoom_scale;

end