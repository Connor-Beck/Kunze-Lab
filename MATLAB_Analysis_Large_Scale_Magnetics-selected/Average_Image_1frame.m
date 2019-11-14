function [AverageImage,ROIbases,Background_Noise] = Average_Image(Image_Stack,num_images,Width,Height,MinVal)

Avg_Image = zeros(Height,Width);


Background_Noise = zeros(num_images,1);

%%
H = waitbar(0,'Creating Average Image');
for k = 1:num_images
    Count = 0;
    waitbar(k/num_images) 
    for j = 1:Width
        for i = 1:Height
            if Avg_Image(i,j) == 0
                if Image_Stack(i,j,k) <= MinVal
                    Avg_Image(i,j) = 0;
                    Background_Noise(k,1) = Image_Stack(i,j,k) + Background_Noise(k,1);
                    Count = Count + 1;
                else
                    Avg_Image(i,j) = 1;
                end
            end
        end
    end
    Background_Noise(k,1) = Background_Noise(k,1)/Count;
end
delete(H)



%%
Image = imbinarize(Image_Stack(:,:,1));
Image = bwareaopen(Image(:,:),40);
AverageImage = imfill(Image(:,:),'holes');

Centroid = regionprops(AverageImage,'centroid');
ROIbases = cat(1, Centroid.Centroid);




end



