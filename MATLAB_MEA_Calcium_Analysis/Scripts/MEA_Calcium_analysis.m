%% MEA-CALCIUM Analysis Script
% Authors: Jack Vincent
%          Andy Kirby
% Last Updated: 2/4/2019 by Andy
%
% Takes the MEA reading and Ca Video as input and detects spikes,
% outputting spike trains as a csv file
%
% Notes:
%   Future plans are to implement NCED detection as opposed to std dev
%   thresholding, as well as restructuring the memory usage of this script

%% Declarations
% MEA sampling rate and paths for channel data and info channel in the MEA 
% hdf file. These may require changing, although they never have for me
Fs = 25000; 
channelds = '/Data/Recording_0/AnalogStream/Stream_0/ChannelData';
labelsds = '/Data/Recording_0/AnalogStream/Stream_0/InfoChannel';

%% Load MEA data and Calcium video
% Open ui to enable user to select hdf and mp4 files and then obtain data 
% for each electrode
[MEA_File, MEA_Path] = uigetfile('../MEAData/*.h5','Select MEA data');
MEA_Data = strcat(MEA_Path, MEA_File); clear MEA_File; clear MEA_Path;
[Calcium_video_file, Calcium_video_path] = uigetfile('../CalciumData/*.mp4','Select Calcium video');
Calcium_video = strcat(Calcium_video_path, Calcium_video_file);
% Get the tiff as well
[Calcium_video_file, Calcium_video_path] = uigetfile('../CalciumData/*.tif','Select Calcium video (tiff)');
Calcium_Video_Filename = strcat(Calcium_video_path, Calcium_video_file);%'2hrnp.tif';%'EcoGlass.tiff';

ChannelData = h5read(MEA_Data,channelds);


%% Preprocess MEA data
% The index of each electrode in the channel data is somewhat random, so 
% this creates a map using the labels field from InfoChannel to convert
% column-row index into linear electrode index in the data
InfoChannel = h5read(MEA_Data,labelsds);
keySet = InfoChannel.Label;
valueSet = 1:60;
labels = containers.Map(keySet,valueSet);

% Obtain electrode of interest column-row index and convert this to the 
% proper linear index in the data using the electrode indexing map
ecrindex = inputdlg('Enter electrode column-row index:');
elindex = labels(ecrindex{1});

% Create time vector for MEA data and extract the data of the EOI
tMEA = 0:1/Fs:(size(ChannelData,1)-1)/Fs;
EOI = ChannelData(:,elindex);

% Obtain conversion factor of data (depends on compression of hdf), then
% apply it and convert the electrode data from picovolts to microvolts
ConversionFactor = double(InfoChannel.ConversionFactor(1));
EOI = 10^-6*ConversionFactor*double(EOI);

%% Load Calcium video, count frames, and establish baseline
% Load calcium video and count the number of frames. Reload the video 
% because MATLAB is weird
v = VideoReader(Calcium_video);
frames = v.NumFrames;
v = VideoReader(Calcium_video);

%Obtain the coordinates of the electrode of interest
firstFrame = readFrame(v);
figure('WindowState','maximized')
imshow(firstFrame);
title('Select Center of Electrode')
[~, xi2, yi2] = roipoly;
close
v.CurrentTime = 0;
electrode_coord = [yi2 xi2];
clear xi2
clear yi2

%%
% Pull up the first frame of the video, and then manually select a 
% background region which will be used to establish a baseline. Reset video
% to t = 0
firstFrame = readFrame(v);
figure('WindowState','maximized')
imshow(firstFrame);
title('Background Region')
baselineRegion = roipoly;
close
v.CurrentTime = 0;

% Initialize baseline vector which will act as F0
baseline = zeros(frames,1);

% Read every frame in the video, take the average of the background region,
% store it in the baseline vector, and then reset the video to t = 0
H = waitbar(0,'Determining baseline fluorescence...');
k = 1;
while hasFrame(v)
    waitbar(k/frames)
    video = readFrame(v);
    video = rgb2gray(video);
    baseline(k) = mean2(video(baselineRegion));
    k = k+1;
end
delete(H)
v.CurrentTime = 0;

%% Detect Calcium ROIs, extract data, and interpolate through anomalies
% Uses @Connor Beck's segmentation methods to determine ROI
%Make sure to use tiff format (can convert in imageJ)
[Image_Stack,num_images,Width,Height] = Image_Reader_editor(Calcium_Video_Filename);

%%
%Adjust MinVal based on base F
%Prompt for MinVal. MinVal corresponds to the minimum "heat" in the average
%image for a region to be counted as a ROI
MinVal = inputdlg('Enter MinVal for F detection (0-255):');
MinVal = str2num(MinVal{1});

%Create an "average image" which essentially finds the "hot spots" within
%the video
AverageImage = Average_Image_editor(Image_Stack,num_images,Width,Height,MinVal);
imshow(AverageImage)
%%

%ROI is a cell array of the yx coordinates of each region, L is something
%(ask Connor, doubt he knows, good thing he comments his code!, ROIboundary
%contains the boundary of the region
[ROI, L, ROIboundary] = segmentation(AverageImage);

%Displays the detected ROI for manual check
%FUTURE: Add way to edit/delete/maybe add ROI manually
ShowROI(Image_Stack(:,:,1),ROIboundary,L)

%%
%---------------Manual ROI Editting----------------------------------------
% ~~~~~~~~~DELETION~~~~~~~~~~~~
% Ask for ROI to delete
while 1==1
    ROI_to_delete = inputdlg('Enter ROI to delete (Cancel to stop deleting):');
    if size(ROI_to_delete) == 0
       break 
    end
    ROI_to_delete = str2num(ROI_to_delete{1});
    ROI(ROI_to_delete) = [];
    ROIboundary(ROI_to_delete) = [];
    
    ShowROI(Image_Stack(:,:,1),ROIboundary,L)
end
%%
% ~~~~~~~~~ADDITION~~~~~~~~~~~~
% Determine if user would like to add more ROIs
answer = questdlg('Would you like to add another ROI?');

% Initalize binary variable to determine if user is done selecting ROIs
ROIsDone = 0;
while ROIsDone == 0
    % As long as the user keeps asking to add more ROIs, pull up the first
    % frame of the video, and add another ROI 
    switch answer
        case 'Yes'
            figure('WindowState','maximized')
            ShowROI(Image_Stack(:,:,1),ROIboundary,L)
            %imshow(firstFrame)
            title('Select a ROI to add')
            tempROImask = roipoly; %Prompt the ROI mask
            [tempROIboundary, tempL] = bwboundaries(tempROImask); %Find the boundary
            %Now find the coordinates of each pixel in the ROI
            pixel = 1;
            tempROI = cell(0);
            for pixC = 1:size(tempROImask,1)
                for pixR = 1:size(tempROImask,2)
                    if tempROImask(pixC,pixR) == 1
                        tempROI{pixel,1} = [pixC,pixR];
                        pixel = pixel+1;
                    end
                end
            end
            ROI{size(ROI,1)} = tempROI;
            ROIboundary{size(ROIboundary,1)} = tempROIboundary{1};
            close
            answer = questdlg('Would you like to add another ROI?');
        case 'No'
            ROIsDone = 1;
        case 'Cancel'
            ROIsDone = 1;
        case ''
            ROIsDone = 1;
    end
end
%clear ROIsDone tempROI tempROImask tempROIboundary tempL pixel pixC pixR

%%
%In order to port functionality from the ROIs detected, create a mask from
%the coordinates in ROI
H = waitbar(0,'Creating ROI Masks...');
numROIs = size(ROI,1);
ROI_MASKS = cell(numROIs, 1);
for rNum = 1:numROIs
    waitbar(rNum/numROIs)
    ROI_mask = zeros(Height, Width);
    %Access using ROI{ROI_NUMBER}{1}(1 for y, 1 for x)
    %Make each coordinate place a 1, leave the rest 0
    for i = 1:size(ROI{rNum}, 1)
        y = ROI{rNum}{i}(1);
        x = ROI{rNum}{i}(2);
        ROI_mask(y, x) = 1;
    end

    %Convert to logical for compatability
    ROI_mask = logical(ROI_mask);
    
    %Store the ROI mask
    %eval(horzcat('ROI', num2str(rNum), ' = ROI_mask;')); %Takes up a lot
    %of workspace!
    ROI_MASKS{rNum, 1} = ROI_mask;
end
delete(H)
clear rNum

%%
% Initialize Calcium time vector and ROI fluorescence vectors to store
% normalized ROI fluorescence for each ROI
tCalcium = 0:1/v.FrameRate:(frames-1)/v.FrameRate;
fluoROI = cell(numROIs, 1);
for i = 1:numROIs
    %eval(horzcat('fluoROI',num2str(i),' = zeros(frames,1);'));
    fluoROI{i,1} = zeros(frames,1);
end

%%
% Determine normalized ROI fluorescence by extracting fluoresence of each 
%  ROI for each frame, subtracting baseline, and then dividing by baseline
%
% Note: Time intensive, maybe look into a method of decreasing this time?
% Blows up with greater number of ROI. Took about 
H = waitbar(0,'Extracting fluorescence of each ROI');
k = 1;
v.CurrentTime = 0;
while hasFrame(v)
    waitbar(k/frames)
    video = readFrame(v);
    video = rgb2gray(video);
    for i = 1:numROIs
        %eval(horzcat('fluoROI',num2str(i),'(k) = (mean2(video(ROI' ...
        %    ,num2str(i),'))-baseline(k))/baseline(k);'));
        fluoROI{i}(k) = (mean2(video(ROI_MASKS{i}))-baseline(k))/baseline(k);
    end
    k = k+1;
end
delete(H)
%%
%Find and interpolate anomalies for every single ROI
j = 1;
while j ~= numROIs + 1
    
    % Initalize binary variable which determines whether or not all 
    % anomalies for this ROI have been eliminated
    anomaliesElim = 0;
    
    % Continue asking to eliminate anomalies until all are eliminated for
    % this ROI
    while anomaliesElim == 0

        % Plot current ROI fluorescence over time
        figure('WindowState','maximized')
        plot(tCalcium, fluoROI{j})%eval(horzcat('fluoROI',num2str(j))))
        xlim([0 v.Duration]);
        title(horzcat('ROI ',num2str(j),' Fluorescence'))
        xlabel('Time (s)')
        ylabel('\DeltaF/F')

        % Determine if user would like to eliminate an anomaly
        answer = questdlg(['Would you like to erase and interpolate ' ...
            'an anomaly?']);

        % If anything other than yes, switch anomalies eliminated variable 
        % to true, so while loop will break on next run. If yes, eliminate
        % the anomaly
        switch answer
            case 'Yes'
                % Give user instructions and initialize binary keyboard 
                % button pressed variable
                title(['Use the zoom function to select the portion ' ...
                    'of the graph you would like to erase and ' ...
                    'interpolate and then press any key to continue'])
                k = 0;
                % As long as keyboard button hasn't been pressed, wait
                while k == 0
                    k = waitforbuttonpress;
                    % When keyboard button is pressed, determine the x 
                    % limits of the current axes, determine the closest 
                    % time indices of these x limits, and create a linear 
                    % function connecting the first point on the left and 
                    % the last point on the right
                    if k == 1
                        ax = gca;
                        deleteTimes = ax.XLim;
                        [~,startIdx] = min(abs(tCalcium-deleteTimes(1)));
                        [~,endIdx] = min(abs(tCalcium-deleteTimes(2)));
                        m = (fluoROI{j}(endIdx)-fluoROI{j}(startIdx))/(tCalcium(endIdx)-tCalcium(startIdx));
                        b = fluoROI{j}(startIdx) - m*tCalcium(startIdx);
                        %eval(horzcat('m = (fluoROI',num2str(j), ...
                        %    '(endIdx)-fluoROI',num2str(j), ...
                        %    '(startIdx))/(tCalcium(endIdx)', ...
                        %    '-tCalcium(startIdx));'));
                        %eval(horzcat('b = fluoROI',num2str(j), ...
                        %   '(startIdx)-m*tCalcium(startIdx);'));
                        
                        % Use this linear function to fill in the gap 
                        % between the lower and uppper x limits of anomaly
                        for i = startIdx:endIdx
                            %eval(horzcat('fluoROI',num2str(j), ...
                            %    '(i) = m*tCalcium(i)+b;'));
                            fluoROI{j}(i) = m*tCalcium(i)+b;
                        end
                        close
                    end 
                end
            case 'No'
                close
                anomaliesElim = 1;
            case 'Cancel' %%Exit from all ROI anomalies
                close
                anomaliesElim = 1;
                j = numROIs;
            case ''
                close
                anomaliesElim = 1;
        end
    end
    j = j + 1;
end

%% MEA spike detection and recording
% Initialize two time vectors: one that will be used as a binary mask to
% denote locations of spikes (tspikesMEA) and one that will be used to
% store spike time stamps for graphing and analysis (xspikesMEA)
tspikesMEA = zeros(size(tMEA,2),1);
xspikesMEA = zeros(size(tMEA,2),1);

% Calculate standard deviation of electrode of interest. Multiply by 5.5 to
% obtain threshold per Ide et al. (2010)
EOIstd = std(EOI);
EOIthresh = -4.75*EOIstd; %-4.5*EOIstd;

% Initialize below threshold binary variable
if EOI(1) < EOIthresh
    belowthresh = 1;
else
    belowthresh = 0;
end

% Comb through EOI data. If value is below threshold and data was
% previously above threshold, mark the time and correct the below threshold 
% variable. If value is above threshold and data was previously below 
% threshold, just correct the below threshold variable.
for i = 1:size(EOI)
    if EOI(i) < EOIthresh && belowthresh == 0
        tspikesMEA(i) = 1;
        belowthresh = 1;
    elseif EOI(i) > EOIthresh && belowthresh == 1
        belowthresh = 0;
    end
end

% Use the binary mask to obtain and store spike times
for i = 1:size(tspikesMEA)
    if tspikesMEA(i) == 1
       xspikesMEA(i) = tMEA(i);
    end
end

% Condense spike times vector by deleting all zeros
xspikesMEA(xspikesMEA == 0) = [];

% Delete all spike timestamps that fall within 2 ms of a previous spike
i = 2;
k = size(xspikesMEA,1);
while i <= k
    if xspikesMEA(i,1) < xspikesMEA(i-1,1) + 0.002
        xspikesMEA(i,:) = [];
    end
    i = i+1;
    k = size(xspikesMEA,1);
end

% Create a second column copy of spike timestamps. This is for graphing
% vertical lines at the spike timestamps with the line() function later
xspikesMEA = horzcat(xspikesMEA,xspikesMEA);

%% Calcium spike detection and recording
% Calculate the first order time derivative of the ROI fluorescence for 
% each ROI. Dividing by seconds per frame is the same thing as multiplying 
% by the frame rate
fluoROIdiff = cell(numROIs, 1);
for i = 1:numROIs
    %eval(horzcat('fluoROIdiff',num2str(i),' = diff(fluoROI',num2str(i), ...
    %    ')*v.FrameRate;'));
    fluoROIdiff{i} = diff(fluoROI{i})*v.FrameRate;
end

% Calculate the new time vector associated with these ROI derivatives
tCalciumStep = tCalcium(2)-tCalcium(1);
tCalciumDiff = tCalciumStep:1/v.FrameRate:tCalcium(end);

% Bin the ROI fluorescence derivative data in groups of 2 for each ROI
for j = 1:numROIs
    %eval(horzcat('fluoROIdiff = fluoROIdiff',num2str(j),';'));
    for i = 1:floor(size(fluoROIdiff{j},1)/2)
        fluoROIdiff{j}(i) = fluoROIdiff{j}(i)+fluoROIdiff{j}(i+1);
        fluoROIdiff{j}(i+1) = [];
    end
    %eval(horzcat('fluoROIdiff',num2str(j),' = fluoROIdiff;'));
end

% Average every 2 timestamps to obtain the corresponding time vector for
% the binned ROIs derivative data. Subtract 500 ms for delay of fluo-4
% binding kinetics
for i = 1:floor(size(tCalciumDiff,2)/2)
    tCalciumDiff(i) = (tCalciumDiff(i)+tCalciumDiff(i+1))/2;
    tCalciumDiff(i+1) = [];
end

tCalciumDiff = tCalciumDiff-0.5;
%%
% Determine spike locations for every ROI
Calcium_Thresholds = zeros(numROIs + 1, 1); % Initialize vector to record the Ca Thresholds
fluoROIdiffstd = cell(numROIs, 1);
xspikesCalciumCells = cell(numROIs, 1);
for j = 1:numROIs
    % Once again, initialize two time vectors: a binary mask and one for
    % calcium spike timestamps 
    tspikesCalcium = zeros(size(tCalciumDiff,2),1);
    xspikesCalcium = zeros(size(tCalciumDiff,2),1);

    % Determine Calcium threshold in accordance with Ide et al.
    %eval(horzcat('fluoROIdiffstd',num2str(j), ...
    %    ' = std(fluoROIdiff',num2str(j),');'));
    %eval(horzcat('CalciumThresh = 2*fluoROIdiffstd',num2str(j)));
    %eval(horzcat('fluoROIdiff = fluoROIdiff',num2str(j),';'));
    fluoROIdiffstd{j} = std(fluoROIdiff{j});
    CalciumThresh = 2*fluoROIdiffstd{j};
    
    %Store the calcium thresholds
    Calcium_Thresholds(j+1, 1) = CalciumThresh;

    % Initalize a binary calcium threshold variable
    if fluoROIdiff{j}(1) < CalciumThresh
        belowthresh = 1;
    else
        belowthresh = 0;
    end

    % Comb through ROI data. If value is above threshold and data was
    % previously below threshold, mark the time and correct the below 
    % threshold variable. If value is below threshold and data was 
    % previously above threshold, just correct the threshold variable.
    for i = 1:size(fluoROIdiff{j})
        if fluoROIdiff{j}(i) > CalciumThresh && belowthresh == 1
            tspikesCalcium(i) = 1;
            belowthresh = 0;
        elseif fluoROIdiff{j}(i) < CalciumThresh && belowthresh == 0
            belowthresh = 1;
        end
    end

    % Use the binary mask to obtain and store spike times
    for i = 1:size(tspikesCalcium)
        if tspikesCalcium(i) == 1
            xspikesCalcium(i) = tCalciumDiff(i);
        end
    end
   
    % Condense spike times vector by deleting all zeros
    xspikesCalcium(xspikesCalcium == 0) = [];

    % Delete all spike timestamps that fall within 5 s of a previous spike
    i = 2;
    k = size(xspikesCalcium,1);
    while i <= k
        if xspikesCalcium(i,1) < xspikesCalcium(i-1,1) + 2
            xspikesCalcium(i,:) = [];
        end
        i = i+1;
        k = size(xspikesCalcium,1);
    end

    % Create a second column copy of spike timestamps. This is for graphing
    % vertical lines at the spike timestamps with the line() function
    % later. Store spike data in the appropriate calcium spikes ROI
    % variable number 
    xspikesCalcium = horzcat(xspikesCalcium,xspikesCalcium);
    %eval(horzcat('xspikesCalcium',num2str(j),' = xspikesCalcium;'));
    xspikesCalciumCells{j} = xspikesCalcium;
end

%% Plot MEA data
% Pull up a new maximized figure
figure('WindowState','maximized')

% Plot MEA data for electrode of interest. Add a red line to mark threshold 
% and give its value in the legend. Add green vertical lines to mark spikes
subplot(2,1,1)
plot(tMEA,EOI)
xlim([0 v.Duration])
hline = refline(0,-4*EOIstd);
hline.Color = 'r';
for i = 1:size(xspikesMEA,1)
    line(xspikesMEA(i,:),get(gca,'YLim'),'Color','g')
end
title(horzcat('Electrode ',ecrindex{1}))
ylabel('Voltage (\muV)')
xlabel('Time (s)')
legend(hline,horzcat('Threshold = ',num2str(EOIthresh),' \muV'))

%% Plot ROI data
% Set ROI you want to look at with the currentROI variable. Set all the 
% necessary plotting variables to the appropriate number of ROI. Beneath 
% the MEA data, plot Calcium first order derivative data. Also add a red 
% line for threshold and give its value in the legend and mark spikes with 
% vertical green lines
currentROI = 1;
%eval(horzcat('fluoROIdiff = fluoROIdiff',num2str(currentROI),';'));
%eval(horzcat('fluoROIdiffstd = fluoROIdiffstd',num2str(currentROI),';'));
%eval(horzcat('xspikesCalcium = xspikesCalcium',num2str(currentROI),';'));

subplot(2,1,2)
plot(tCalciumDiff,fluoROIdiff{currentROI})
xlim([0 v.Duration])
hline = refline(0,fluoROIdiffstd{currentROI});
% = refline(0,num2str(fluoROIdiffstd{currentROI}));
hline.Color = 'r';
for i = 1:size(xspikesCalciumCells{currentROI},1)
    %line(xspikesCalcium{currentROI}(i,:),get(gca,'YLim'),'Color','g')
end
%title(horzcat('ROI ',num2str(currentROI), ...
    %' Fluorescence (1st Order Time Derivative)'))
title(horzcat('ROI ',num2str(currentROI), ...
    ' Fluorescence'))
xlabel('Time (s)')
ylabel('d(\DeltaF/F)/dt')
legend(hline,horzcat('Threshold = ',num2str(fluoROIdiffstd{currentROI})))

%% Outputs
% Ouput the number of concurrent spike events, their timestamps, the ROIs 
% associated with the Calcium spikes, the average electrode firing rate and 
% the number of noncurrent MEA and Calcium spike events

disp(['Electrode firing rate: ' ... 
    num2str(60*(size(xspikesMEA,1)/v.Duration)) ' spikes per minute'])

%% Output .csv files for more analysis
% Note: current functionality does not have output table titles,
% should add in later.
%
%Output timestamps of MEA spikes
dlmwrite("MEA_Spikes_" + ecrindex + ".csv", xspikesMEA(:,1), 'precision','%.8f');
%Out timestamps of ROI Calcium Spikes
%Formatted as RO1 RO2 RO3 etc columns

%Create a matrix the size of the max num of spikes and the num of ROIs
caSize = 0;
for i = 1:numROIs-1
    %eval(horzcat('xspikesCalcium = xspikesCalcium',num2str(i),';'));
    if size(xspikesCalciumCells{i}, 1) > caSize
        caSize = size(xspikesCalciumCells{i}, 1);
    end
end
calciumSpikes = zeros(caSize,numROIs);
%Fill the matrix with appropriate values. Will trail with 0s!!
for i = 1:numROIs
    %eval(horzcat('xspikesCalcium = xspikesCalcium',num2str(i),';'));
    for j = 1:size(xspikesCalciumCells{i})
        if size(xspikesCalciumCells{i},2) > 0
            calciumSpikes(j,i) = xspikesCalciumCells{i}(j); 
        end
    end
end
%Write the matrix of Calcium spikes to a file
%Start with the header
fid = fopen('Calcium_Spikes.csv','w');
for i = 1:numROIs
    header={"ROI " + num2str(i)};
    if i == 1
        fprintf(fid,'%s',header{:});
    else
        fprintf(fid,',%s',header{:});
    end
end
fprintf(fid,'\n');
fclose(fid);
%Write data
dlmwrite("Calcium_Spikes.csv", calciumSpikes, '-append');
%---------------------------------------------------------------
%Create table containing number of spikes per bin
bin_size = 2; %seconds
num_bins = size(tCalciumDiff, 2);
bin_data = zeros(num_bins, 2+numROIs);
bin_matrix_header = cell(1, 2+numROIs); 
bin_matrix_header(1) = cellstr("Time");
bin_matrix_header(2) = cellstr("MEA");
for i = 1:numROIs
   bin_matrix_header(i+2) = cellstr("ROI " + num2str(i));
end
bin_matrix_header = cell2table(bin_matrix_header);

cur_col = 1; %keep track of current column we are editing
%Create time bin column
for i = 1:num_bins
   bin_data(i, cur_col) = bin_size * (i-1); 
end
%Increment MEA column (count all fires and place in bin)
%Start with first spike timestamp and first bin, incrementing bins
%once the timestamp passes the bin time, and spikes after they are
%counted
cur_col = cur_col + 1;
cur_bin = 1;
i = 1;
while i <= size(xspikesMEA, 1)
    bin_val = cur_bin * bin_size; %value of time bin we are concerned with
    if xspikesMEA(i, 1) < bin_val %Increment count in the current bin if the timestamp
                                  %lies in it
        bin_data(cur_bin, cur_col) = bin_data(cur_bin, cur_col) + 1;
        i = i + 1;
    else %value must be in next bin
        cur_bin = cur_bin + 1;
        %Make sure we retry counting the current timestamp
    end
end

%Increment ROI columns (count all fires and place in bin)
%Start with first spike timestamp and first bin, incrementing bins
%once the timestamp passes the bin time, and spikes after they are
%counted
for c = 1:numROIs
    cur_col = cur_col + 1;
    cur_bin = 1;
    i = 1;
    while i <= size(calciumSpikes, 1)
        %Break if 0s have been reached
        if calciumSpikes(i, c) == 0 
            break; 
        end
        bin_val = cur_bin * bin_size; %value of time bin we are concerned with
        if calciumSpikes(i, c) < bin_val %Increment count in the current bin if the timestamp
                                      %lies in it
            bin_data(cur_bin, cur_col) = bin_data(cur_bin, cur_col) + 1;
            i = i + 1;
        else %value must be in next bin
            cur_bin = cur_bin + 1;
            %Make sure we retry counting the current timestamp
        end
    end
end
%Write the data table to a file
%Start with the header
fid = fopen('binned_spikes.csv','w');
header={'TIME','MEA'};
fprintf(fid,'%s,%s',header{:});
for i = 1:numROIs
    header={"ROI " + num2str(i)};
    fprintf(fid,',%s',header{:});
end
fprintf(fid,'\n');
fclose(fid);
%Write the data
dlmwrite("binned_spikes.csv", bin_data, '-append');

%Write the EOI Threshold and Ca Thresholds to file (just in case)
Calcium_Thresholds(1) = EOIthresh; %plop that in there to save 
dlmwrite("EOI_ROI_Thresholds.csv", Calcium_Thresholds);

clear header
clear fid
clear bin_matrix_header

%%
%Clear useless variables. Bruh on the legacy code..
clear i
clear j
clear k
clear ans
clear answer
clear c
clear CaSize
clear H
clear r
clear x
clear y
clear belowthresh






