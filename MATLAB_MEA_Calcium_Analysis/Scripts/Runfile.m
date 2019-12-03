

[Image_Stack,num_images,Width,Height] = Image_Reader_editor('EcoGlass.tiff');
%Use tif format (Can convert in imageJ)

%% Locate ROIs for Analysis


MinVal =  180;
%Need to set some sort of min val


AverageImage = Average_Image_editor(Image_Stack,num_images,Width,Height,MinVal);
[ROI, L, ROIboundary] = segmentation(AverageImage);
%ROI cell array with the xy coordinates of each region

%%
%Check to make sure ROIs look good

ShowROI(Image_Stack(:,:,1),ROIboundary,L)


%%
%Create mask from the coordinates
ROI_MASK = zeros(Height, Width);
%Access using ROI{ROI_NUMBER}{1}(1 for y, 1 for x)
num_coord = size(ROI{1}, 1);

for i = 1:num_coord
    y = ROI{1}{i}(1);
    x = ROI{1}{}(2);
    ROI_MASK(y, x) = 1;
end

%Convert to logical
ROI_MASK = logical(ROI_MASK);



%%
%
[DeltaFoverF] = location_analysis_dFoverF(ROI, Image_Stack,num_images);


%
Spikes = Spike_Detector(DeltaFoverF);

Show_Spikes(Spikes)


