function [Image_StackTL,Image_StackTR,Image_StackBL,Image_StackBR,num_images,Width,Height] = Image_Reader(filename)

%Obtain values for the image stack that will be used throughout the code
%Save the values in a mat file 'Parameters' for easy loading

ImageNorm = imread(filename);
tiff_info = imfinfo(filename);
num_images = numel(tiff_info);


Width = tiff_info.Width;
Height = tiff_info.Height;


%%Load the filestack and store it in a matrix format
Width = 1024;
Height = 1024;

Image_StackTL = zeros(Height,Width,num_images);
Image_StackTR = zeros(Height,Width,num_images);
Image_StackBL = zeros(Height,Width,num_images);
Image_StackBR = zeros(Height,Width,num_images);

H = waitbar(0,'Loading Images');
for k = 1:num_images
    waitbar(k/num_images)
    Image = imread(filename,'Index', k);
    IMGTL = Image(1:1024,1:1024,1);
    IMGTR = Image(1:1024,1025:2048,1);
    IMGBL = Image(1025:2048,1:1024,1);
    IMGBR = Image(1025:2048,1025:2048,1);
    Image_StackTL(:,:,k) = IMGTL;
    Image_StackTR(:,:,k) = IMGTR;
    Image_StackBL(:,:,k) = IMGBL;
    Image_StackBR(:,:,k) = IMGBR;
end
delete(H)

