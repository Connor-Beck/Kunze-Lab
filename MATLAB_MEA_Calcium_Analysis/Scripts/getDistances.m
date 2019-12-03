%%
%Outputs an array with the distances, in um, that the regions are from the
%point provided. Currently, just uses a random point on the ROI boundary.
%Requires that the scale is given (in um per pixel).
function [DISTANCES] = getDistances(point, ROIboundary, um_per_pixel)
    DISTANCES = double(size(ROIboundary,1));
    p_x = point(1);
    p_y = point(2);
    for i = 1:size(ROIboundary,1)
       x = ROIboundary{i}(1,1);
       y = ROIboundary{i}(1,2);
       diff_x = p_x - x;
       diff_y = p_y - y;
       pixel_dist = sqrt(diff_x^2 + diff_y^2);
       DISTANCES(i) = pixel_dist * um_per_pixel;
    end
end