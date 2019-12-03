%% Declarations
% MEA sampling rate and paths for channel data and info channel in the MEA 
% hdf file. These may require changing, although they never have for me
Fs = 25000; 
channelds = '/Data/Recording_0/AnalogStream/Stream_0/ChannelData';
labelsds = '/Data/Recording_0/AnalogStream/Stream_0/InfoChannel';

%% Load MEA data and Calcium video
% Open ui to enable user to select hdf and mp4 files and then obtain data 
% for each electrode
MEA_Data = uigetfile('*.h5','Select MEA data');
Calcium_video = uigetfile('.mp4','Select Calcium video');
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
frames = v.NumberofFrames;
v = VideoReader(Calcium_video);

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
k = 1;
while hasFrame(v)
    video = readFrame(v);
    video = rgb2gray(video);
    baseline(k) = mean2(video(baselineRegion));
    k = k+1;
end
v.CurrentTime = 0;

%% Select Calcium ROIs, extract data, and interpolate through anomalies

% Export and open first frame of video in default external software, so 
% user can use Paint or some other software to plan out their ROIs, wait 
% until they do anything with the dialog box to continue
firstFramefilename = [Calcium_video ' (FIRST FRAME).jpg'];
imwrite(firstFrame,firstFramefilename)
answer = questdlg(['First frame of video exported for ROI planning.' ...
    'Press any button to continue.']);
switch answer
    case 'Yes'
    case 'No'
    case 'Cancel'
    case ''
end

% Pull up the first frame of the video again and select the neuron whose 
% fluorescence you're interested in tracking
figure('WindowState','maximized')
imshow(firstFrame)
title('ROI 1')
ROI1 = roipoly;
close

% Determine if user would like to add more ROIs
numROIs = 1;
answer = questdlg('Would you like to add another ROI?');

% Initalize binary variable to determine if user is done selecting ROIs
ROIsDone = 0;
while ROIsDone == 0
    % As long as the user keeps asking to add more ROIs, pull up the first
    % frame of the video, and add another ROI 
    switch answer
        case 'Yes'
            figure('WindowState','maximized')
            imshow(firstFrame)
            title(horzcat('ROI ',num2str(numROIs+1)))
            tempROI = roipoly;
            eval(horzcat('ROI',num2str(numROIs+1),' = tempROI;'));
            close
            numROIs = numROIs+1;
            answer = questdlg('Would you like to add another ROI?');
        case 'No'
            ROIsDone = 1;
        case 'Cancel'
            ROIsDone = 1;
        case ''
            ROIsDone = 1;
    end
end

% Initialize Calcium time vector and ROI fluorescence vectors to store
% normalized ROI fluorescence for each ROI
tCalcium = 0:1/v.FrameRate:(frames-1)/v.FrameRate;
for i = 1:numROIs
    eval(horzcat('fluoROI',num2str(i),' = zeros(frames,1);'));
end

% Determine normalized ROI fluorescence by extracting fluoresence of each 
%  ROI for each frame, subtracting baseline, and then dividing by baseline
k = 1;
while hasFrame(v)
    video = readFrame(v);
    video = rgb2gray(video);
    for i = 1:numROIs
        eval(horzcat('fluoROI',num2str(i),'(k) = (mean2(video(ROI' ...
            ,num2str(i),'))-baseline(k))/baseline(k);'));
    end
    k = k+1;
end

%Find and interpolate anomalies for every single ROI
for j = 1:numROIs
    
    % Initalize binary variable which determines whether or not all 
    % anomalies for this ROI have been eliminated
    anomaliesElim = 0;
    
    % Continue asking to eliminate anomalies until all are eliminated for
    % this ROI
    while anomaliesElim == 0

        % Plot current ROI fluorescence over time
        figure('WindowState','maximized')
        plot(tCalcium,eval(horzcat('fluoROI',num2str(j))))
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
                        eval(horzcat('m = (fluoROI',num2str(j), ...
                            '(endIdx)-fluoROI',num2str(j), ...
                            '(startIdx))/(tCalcium(endIdx)', ...
                            '-tCalcium(startIdx));'));
                        eval(horzcat('b = fluoROI',num2str(j), ...
                            '(startIdx)-m*tCalcium(startIdx);'));
                        % Use this linear function to fill in the gap 
                        % between the lower and uppper x limits of anomaly
                        for i = startIdx:endIdx
                            eval(horzcat('fluoROI',num2str(j), ...
                                '(i) = m*tCalcium(i)+b;'));
                        end
                        close
                    end 
                end
            case 'No'
                close
                anomaliesElim = 1;
            case 'Cancel'
                close
                anomaliesElim = 1;
            case ''
                close
                anomaliesElim = 1;
        end
    end
end

%% MEA spike detection and recording
% Initialize two time vectors: one that will be used as a binary mask to
% denote locations of spikes (tspikesMEA) and one that will be used to
% store spike time stamps for graphing and analysis (xspikesMEA)
tspikesMEA = zeros(size(tMEA,2),1);
xspikesMEA = zeros(size(tMEA,2),1);

[xspikesMEA, tspikesMEA, EOItest, NCED] = ncedDetect(EOI, xspikesMEA, tspikesMEA, Fs);
xspikesMEA = horzcat(xspikesMEA,xspikesMEA);

%% Calcium spike detection and recording
% Calculate the first order time derivative of the ROI fluorescence for 
% each ROI. Dividing by seconds per frame is the same thing as multiplying 
% by the frame rate
for i = 1:numROIs
    eval(horzcat('fluoROIdiff',num2str(i),' = diff(fluoROI',num2str(i), ...
        ')*v.FrameRate;'));
end

% Calculate the new time vector associated with these ROI derivatives
tCalciumStep = tCalcium(2)-tCalcium(1);
tCalciumDiff = tCalciumStep:1/v.FrameRate:tCalcium(end);

% Bin the ROI fluorescence derivative data in groups of 2 for each ROI
for j = 1:numROIs
    eval(horzcat('fluoROIdiff = fluoROIdiff',num2str(j),';'));
    for i = 1:floor(size(fluoROIdiff,1)/2)
        fluoROIdiff(i) = fluoROIdiff(i)+fluoROIdiff(i+1);
        fluoROIdiff(i+1) = [];
    end
    eval(horzcat('fluoROIdiff',num2str(j),' = fluoROIdiff;'));
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
for j = 1:numROIs
    % Once again, initialize two time vectors: a binary mask and one for
    % calcium spike timestamps 
    tspikesCalcium = zeros(size(tCalciumDiff,2),1);
    xspikesCalcium = zeros(size(tCalciumDiff,2),1);

    % Determine Calcium threshold in accordance with Ide et al.
    eval(horzcat('fluoROIdiffstd',num2str(j), ...
        ' = std(fluoROIdiff',num2str(j),');'));
    eval(horzcat('CalciumThresh = 2*fluoROIdiffstd',num2str(j)));
    eval(horzcat('fluoROIdiff = fluoROIdiff',num2str(j),';'));
    
    %Store the calcium thresholds
    Calcium_Thresholds(j+1, 1) = CalciumThresh;

    % Initalize a binary calcium threshold variable
    if fluoROIdiff(1) < CalciumThresh
        belowthresh = 1;
    else
        belowthresh = 0;
    end

    % Comb through ROI data. If value is above threshold and data was
    % previously below threshold, mark the time and correct the below 
    % threshold variable. If value is below threshold and data was 
    % previously above threshold, just correct the threshold variable.
    for i = 1:size(fluoROIdiff)
        if fluoROIdiff(i) > CalciumThresh && belowthresh == 1
            tspikesCalcium(i) = 1;
            belowthresh = 0;
        elseif fluoROIdiff(i) < CalciumThresh && belowthresh == 0
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
    eval(horzcat('xspikesCalcium',num2str(j),' = xspikesCalcium;'));
end

%% Plot MEA data
% Pull up a new maximized figure
figure('WindowState','maximized')

% Plot MEA data for electrode of interest. Add a red line to mark threshold 
% and give its value in the legend. Add green vertical lines to mark spikes
subplot(2,1,1)
plot(tMEA,EOI)
xlim([0 v.Duration])
%hline = refline(0,-4*EOIstd);
%hline.Color = 'r';
for i = 1:size(xspikesMEA,1)
    %line(xspikesMEA(i,:),get(gca,'YLim'),'Color','g')
end
title(horzcat('Electrode ',ecrindex{1}))
ylabel('Voltage (\muV)')
xlabel('Time (s)')
%legend(hline,horzcat('Threshold = ',num2str(EOIthresh),' \muV'))

%% Plot ROI data
% Set ROI you want to look at with the currentROI variable. Set all the 
% necessary plotting variables to the appropriate number of ROI. Beneath 
% the MEA data, plot Calcium first order derivative data. Also add a red 
% line for threshold and give its value in the legend and mark spikes with 
% vertical green lines
currentROI = 1;
eval(horzcat('fluoROIdiff = fluoROIdiff',num2str(currentROI),';'));
eval(horzcat('fluoROIdiffstd = fluoROIdiffstd',num2str(currentROI),';'));
eval(horzcat('xspikesCalcium = xspikesCalcium',num2str(currentROI),';'));

subplot(2,1,2)
plot(tCalciumDiff,fluoROIdiff)
xlim([0 v.Duration])
%--hline = refline(0,fluoROIdiffstd);
%hline = refline(0,num2str(fluoROIdiffstd));
%hline.Color = 'r';
for i = 1:size(xspikesCalcium,1)
    %line(xspikesCalcium(i,:),get(gca,'YLim'),'Color','g')
end
%title(horzcat('ROI ',num2str(currentROI), ...
    %' Fluorescence (1st Order Time Derivative)'))
title(horzcat('ROI ',num2str(currentROI), ...
    ' Fluorescence'))
xlabel('Time (s)')
ylabel('d(\DeltaF/F)/dt')
%legend(hline,horzcat('Threshold = ',num2str(fluoROIdiffstd)))

%% Check for concurrent spikes
% Initalize a structure with a field to store concurrent MEA spikes, a 
% field to store concurrent Calcium spike timestamps, and a field to store
% the ROI of every concurrent Calcium spike
concurrentSpikes = cell(size(tCalciumDiff,2),4);
for i = 1:size(tCalciumDiff,2)
    concurrentSpikes{i,1} = tCalciumDiff(i);
end

for i = 1:size(xspikesMEA,1)
    [~,idx] = min(abs(tCalciumDiff - xspikesMEA(i)));
    concurrentSpikes{idx,2} = [concurrentSpikes{idx,2} xspikesMEA(i)];
end

for i = 1:numROIs
    eval(['xspikesCalcium = xspikesCalcium' num2str(i) ';']);
    for j = 1:size(xspikesCalcium,1)
        [~,idx] = min(abs(tCalciumDiff - xspikesCalcium(j)));
        concurrentSpikes{idx,3} = [concurrentSpikes{idx,3} ...
            xspikesCalcium(j)];
        concurrentSpikes{idx,4} = [concurrentSpikes{idx,4} i];
    end
end

numCalcium = 0;
for i = 1:numROIs
    eval(['xspikesCalcium = xspikesCalcium' num2str(i) ';']);
    numCalcium = numCalcium + size(xspikesCalcium,1);
end

numNC = 0;
for i = 1:size(tCalciumDiff,2)  
    numNCthisTime = (size(concurrentSpikes{i,2},2) ...
        - size(concurrentSpikes{i,3},2));
    if  numNCthisTime < 0
        numNCthisTime = 0;
    end
    numNC = numNC + numNCthisTime;
end

numMEA_NC = numNC;

numNC = 0;
for i = 1:size(tCalciumDiff,2)
    numNCthisTime = (size(concurrentSpikes{i,3},2) ...
        - size(concurrentSpikes{i,2},2));
    if  numNCthisTime < 0
        numNCthisTime = 0;
    end
    numNC = numNC + numNCthisTime;
end

numCalcium_NC = numNC;

%% Outputs
% Ouput the number of concurrent spike events, their timestamps, the ROIs 
% associated with the Calcium spikes, the average electrode firing rate and 
% the number of noncurrent MEA and Calcium spike events

disp(['Electrode firing rate: ' ... 
    num2str(60*(size(xspikesMEA,1)/v.Duration)) ' spikes per minute'])

disp(horzcat('Number of nonconcurrent MEA spike events: ', ...
    num2str(numMEA_NC)))

disp(horzcat('Number of nonconcurrent Calcium spike events: ', ...
    num2str(numCalcium_NC)))

disp(['Percentage of Calcium spikes accounted for: ' ...
    num2str(100-(numCalcium_NC/numCalcium)*100) '%'])

numEmpty = 0;
numFull = 0;
for i = 1:size(concurrentSpikes,1)
    if isempty(concurrentSpikes{i,2})
        numEmpty = numEmpty+1;
    else
        numFull = numFull+1;
    end
end

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
    eval(horzcat('xspikesCalcium = xspikesCalcium',num2str(i),';'));
    if size(xspikesCalcium, 1) > caSize
        caSize = size(xspikesCalcium, 1);
    end
end
calciumSpikes = zeros(caSize,numROIs);
%Fill the matrix with appropriate values. Will trail with 0s!!
for i = 1:numROIs
    eval(horzcat('xspikesCalcium = xspikesCalcium',num2str(i),';'));
    for j = 1:size(xspikesCalcium)
        calciumSpikes(j,i) = xspikesCalcium(j); 
    end
end
%Write the matrix of Calcium spikes to a file
dlmwrite("Calcium_Spikes.csv", calciumSpikes);
%---------------------------------------------------------------
%Create table containing number of spikes per bin
bin_size = 2; %seconds
num_bins = size(tCalciumDiff, 2);
bin_data = zeros(num_bins, 2+numROIs);
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
dlmwrite("binned_spikes.csv", bin_data);

%Write the EOI Threshold and Ca Thresholds to file (just in case)
Calcium_Thresholds(1) = EOIthresh; %plop that in there to save 
dlmwrite("EOI_ROI_Thresholds.csv", Calcium_Thresholds);




