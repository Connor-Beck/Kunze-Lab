function Optimal_MinVal = Optimize_MinVal(Image_Stack,num_images,Width,Height)

maxROI=zeros(1);
i=60;
while i < 150
    MinVal = i;
    [AverageImage,ROIbases,BackgroundNoise] = Average_Image(Image_Stack,num_images,Width,Height,MinVal);
    [ROI, L, ROIboundary] = segmentationNEW(AverageImage);   
    if length(ROI) >= length(maxROI)
        maxAverageImage=AverageImage;
        maxROIbases=ROIbases;
        maxBackgroundNoise=BackgroundNoise;
        maxROI=ROI;
        maxL=L;
        maxROIboundary=ROIboundary;
        Optimal_MinVal=i;
    end
    i=i+10;
end

ShowROIedits(maxAverageImage,maxROIboundary,maxL)