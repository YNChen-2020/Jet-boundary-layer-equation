function [centerline_index,edge_index] = find_Edge(velocity_magnitude,zoom_scale,ratio)

% find the max velocity of the velocity data and the edge of the velocity
% the velolocity_magnitude is the velocity data
% the zoom_scale is the Interpolation expansion factor
% the ratio is the edge_velocity to the max_velocity

% delete the small and wrong data
velocity_result = velocity_magnitude;
velocity_result(velocity_result<1e-3) = 0;       %low pass the small value
find_edge = edge(velocity_result,'Canny');       %find the edge and the range of jets
find_null = sum(find_edge) == 0;
velocity_result(:, find_null) = [];

velocity_Zoom = imresize(velocity_result, zoom_scale);
[row, column] = size(velocity_Zoom);
[max_index, edge_index1, edge_index2] = deal(zeros(column, 1));
for j = 1 : column
    velocity_column = velocity_Zoom(:, j);
    max_velocity = max(velocity_column);
    edge_velocity = ratio * max_velocity;         % denote the value of velocity on the edge
    velocity_diff = abs(velocity_column - edge_velocity);
    max_index(j, 1) = find(velocity_column == max_velocity);
    [~, edge_index1(j,1)] = min(velocity_diff(1 : max_index(j, 1)));         % find the edge of r (1/2)
    [~, edge_index2(j,1)] = min(velocity_diff(max_index(j, 1) : end));
    edge_index2(j, 1) = edge_index2(j, 1) + max_index(j, 1) - 1;
end

% zoom out the index to the original data (size and the value)
max_index_zoomout = imresize(velocity_result, 1/zoom_scale);
centerline_index = (row - max_index_zoomout) ./ zoom_scale;
edge_inde_zoom_out = [imresize(edge_index1, 1/zoom_scale),...
    imresize(edge_index2, 1/zoom_scale)];
edge_index = (row - edge_inde_zoom_out) ./ zoom_scale;

end




