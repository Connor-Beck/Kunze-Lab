clear
clc
%function [dF,Spikes,ROI,ROIbases,corr] =
%ConnectivityAnalysis(filename,MinValROI,std_threshold,static_threshold,CellMapFig)

% filename='Videos/Compressed/04_wNP_wmag.tiff';
% MinValROI = 1000;
% std_threshold = 10;
% static_threshold = 0.01;
% CellMapFig=1;

[Image_Stack,num_images,Width,Height] = Image_Reader(filename);

%%
[AverageImage,ROIbases,BackgroundNoise] = Average_Image(Image_Stack,num_images,Width,Height,MinValROI);
[ROI, L, ROIboundary] = segmentationNEW(AverageImage);
%figure(1); ShowROIedits(AverageImage,ROIboundary,L)

dF = dFdtnBKG(ROI, Image_Stack,num_images);
%figure(2); DeltaFoverFplotter(dF,std_threshold,static_threshold)

Spikes = Spike_Detector(dF,std_threshold,static_threshold);
%Spikes(:,3) = 0;

figure(40); spikeplotter = Show_Spikes(Spikes);
close Figure 40

Spike_Count=sum(Spikes);
%Spike_Count=zeros(1,347); Spike_Count(1,1)=10;
Spike_Rate = (Spike_Count/size(Spikes,1));
Mean_Spike_Rate = mean(Spike_Rate)/0.4*60; %spikes/frame*frame/second*second/min 

%close all
TotalROI=size(ROI,1);
%%
corr = correlation_dice(Spikes);
Connected_ROI = Connectivity_dice(corr, ROI);
figure(CellMapFig); Cell_Map_Dice(AverageImage,Connected_ROI,ROIbases);


