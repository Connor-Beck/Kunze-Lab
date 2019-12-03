%displays the ROI 
function showROIbyRadius(AverageImage,ROI,L, ROI_distances, radius)
RASTERFIG = figure('WindowState','maximized');
imshow(rgb2gray(AverageImage));%mat2gray(AverageImage)); 
%imshow(AverageImage);
hold on;
colors=['b' 'g' 'r' 'c' 'm' 'y'];
for k=1:length(ROI)
    if ROI_distances(k) <= radius
        boundary = ROI{k};
        %cidx = mod(color,length(colors))+1;
        cidx = mod(floor(ROI_distances(k)/20),length(colors))+1;
        if(ROI_distances(k) < 140)
            cidx = 1;
        end
        plot(boundary(:,2), boundary(:,1),...
        colors(cidx),'LineWidth',2);
        rndRow = ceil(length(boundary)/(mod(rand*k,7)+1));
        col = boundary(rndRow,2); row = boundary(rndRow,1);
        h = text(col+1, row-1, num2str(L(row,col)));
        set(h,'Color',colors(cidx),'FontSize',14,'FontWeight','bold');
    end
end