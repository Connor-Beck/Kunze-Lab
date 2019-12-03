function [ROI, L, ROIboundary] = segmentation(AverageImage)
% 
% Width = size(AverageImage,2);
% Height = size(AverageImage,1);

[ROIboundary,L,~,~] = bwboundaries(AverageImage);


SingleBoundary = zeros(size(AverageImage,1),size(AverageImage,2),length(ROIboundary));
Cell = 1;
H = waitbar(0,'Segmenting Regions of Interest');

FilledCell = zeros(size(AverageImage,1),size(AverageImage,2),length(ROIboundary));
for i = 1:length(ROIboundary)
    waitbar(i/length(ROIboundary)) 
    Region = ROIboundary{i};
    pixel = 1;
    for j = 1:length(Region)
        SingleBoundary(Region(j,1),Region(j,2),i) = 1;
    end
    FilledCell(:,:,i) = imfill(SingleBoundary(:,:,i),'holes');
    for m = 1:size(AverageImage,1)
        for n = 1:size(AverageImage,2)
            if FilledCell(m,n,i) == 1
                Alpha{pixel,1} = [m,n];
%                 ROI{Cell,1} = Alpha{pixel,1};
                pixel = pixel+1;
            end
        end
    end
    
    ROI{Cell,1} = Alpha;
    Alpha={};
    Cell = Cell+1;
end

delete(H)

% count = 0;
% for i = 1:Height
%     for j = 1:Width
%         if AverageImage(i,j) ~= 0
%             for k = 1:length(ROI)
%                 if i == ROI{:,1} && j == ROI{:,2}
%                     count = count+1;
%                 end
%             end
%         end
%     end
% end
% 
