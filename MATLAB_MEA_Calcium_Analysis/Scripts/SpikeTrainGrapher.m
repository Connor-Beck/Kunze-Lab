%% Spike Train Grapher
%@author Andy Kirby
%Takes input of bin_data and presents a graphical spike train
%representation


%%-------------PLOT MEA AND ALL ROI----------------------------------------
duration = 300;
output_ROI_train_filename = horzcat('spiketrainROIs.jpg');

clear trains
trains{1} = [0];%xspikesMEA(:,[1]);
for i = 1:size(xspikesCalciumCells, 1)
    if size(xspikesCalciumCells{i}, 2) > 0
        trains{i+1} = xspikesCalciumCells{i}(:,[1]);
    end
end

clear labels
labels(1) = "";
%labels(1) = "MEA";
for i = 1:size(xspikesCalciumCells, 1)
    labels(i+1) = "";%horzcat("ROI " + num2str(i)); 
end
%labels(floor(i/2)) = "ROI";
graphTrains(trains, duration, labels, 5);
xlabel("Time (sec)");
xlim([0 200]);
%set(gca,'YColor','RED') % hide the y axis
%ylabel('MEA and ROI')

saveas(gcf, output_ROI_train_filename);

%%
%--------------RADIUS-BASED SPIKETRAIN PLOTTING-------------------------------
%possibly include import ability
%input or import ROI distances/radii
ROI_distances = getDistances(electrode_coord, ROIboundary, 0.426727); %[62,54,70,54,58,62,70,87,74,82,111,120,136,144,185,243,268,338,309];

start_rad = 60;
max_rad = 600;
radii_interval = 10;

num_radii = (max_rad - start_rad) / radii_interval;
radii_to_test = zeros(1, num_radii);
for i = 1:num_radii
    radii_to_test(1, i) = start_rad + (radii_interval*(i-1));
end

radii_to_test = [130, 200, 400];%[130];
%radii_to_test = [85, 90, 100, 105, 110, 115, 120, 180];%[140, 150, 160, 170, 180, 190, 200];

clear trains
trains{1} = xspikesMEA(:,[1]);
%trains{1} = xspikesMEA54_3(:,[1]);
for k = 1:size(radii_to_test, 2)

    cur_radius = radii_to_test(1, k); %in um

    %RASTERFIG2 = figure('WindowState','maximized');
    %showROIbyRadius(Image_Stack(:,:,1),ROIboundary,L, ROI_distances, cur_radius);
    
    %Create the array to hold the Ca spikes
    xspikesCaRadius = zeros(size(tCalciumDiff,2),1);
    indexCaSpikesRadius = 1;
    %Load xspikesCaRadius with appropriate spikes
    for i = 1:size(xspikesCalciumCells, 1)
        if ROI_distances(i) <= cur_radius
            for j = 1:size(xspikesCalciumCells{i}, 1)%size(eval(horzcat('xspikesCalcium', num2str(i))),1)
                if size(xspikesCalciumCells{i},2) > 0
                    %xspikesCaRadius(indexCaSpikesRadius) = eval(horzcat('xspikesCalcium', num2str(i),'(',num2str(j),');'));
                    xspikesCaRadius(indexCaSpikesRadius) = xspikesCalciumCells{i}(j);
                    indexCaSpikesRadius = indexCaSpikesRadius + 1;
                end
            end
        end
    end

    %Clear out empty spots
    xspikesCaRadius(xspikesCaRadius == 0) = [];

    trains{k+1} = xspikesCaRadius;
end
    trains{k+2} = xspikesMEA(:,[1]);
    clear labels
    %labels(1) = "";
    %labels(2) = "";
    labels(1) = "Electrode 54";
    for i = 1:size(radii_to_test,2)
       labels(i+1) = "" + num2str(radii_to_test(1, i)) + "\mum radius"; 
    end
    labels(i+2) = "Electrode 55";
    
    graphTrains(trains, duration, labels, 20);
    %ylabel("ROI in radius 130\mum           MEA");
    xlim([0 200]);
    xlabel("Time (sec)");
    title("Electrode vs Calcium Spiketrains");
    output_RADII_train_filename = horzcat('spiketrain_radii.jpg');
    saveas(gcf, output_RADII_train_filename);
%%
%    RASTERFIG2 = figure('WindowState','maximized');
%    showROIbyRadius(Image_Stack(:,:,1),ROIboundary,L, ROI_distances, 80);%cur_radius);
    % temptrain = trains{15};
    % clear trains
    % trains{1} = xspikesMEA(:,[1]);
    % trains{2} = temptrain
    % graphTrains(trains, duration, ["MEA", "Radius 130"], 40);
    %imshow(firstFrame)
   % showROIbyRadius(firstFrame,ROIboundary,L,ROI_distances,20)
   
   %%
   densehist = zeros(300,1);
   temphist = zeros(300,1);
   temphist2 = zeros(50,1);
   for i = 1:size(xspikesMEA)
      index = floor (xspikesMEA(i));
      temphist(index) = temphist(index) + 1;
   end
   
   for i = 1:50
       for j = 1:6
            index = 6*(i-1) + j;
            temphist2(i) = temphist2(i) + temphist(index);
       end
   end
   
   for i = 1:50
    densehist(i*6) = temphist2(i);
   end
   
   figure
   subplot(2,1,1)
   %plot(temphist2)
   plot(densehist, 'LineWidth', 1.5, 'Color', 'Black')
   ylim([0 10]);
   xlabel('Time (sec)');
   ylabel('Number of spikes');
   title('Electrode 54 Spikes');
   
   densehist = zeros(300,1);
   temphist = zeros(300,1);
   temphist2 = zeros(50,1);
   for i = 1:size(xspikesCaRadius)
      index = floor (xspikesCaRadius(i));
      temphist(index) = temphist(index) + 1;
   end
   
   for i = 1:50
       for j = 1:6
            index = 6*(i-1) + j;
            temphist2(i) = temphist2(i) + temphist(index);
       end
   end
   
   for i = 1:50
    densehist(i*6) = temphist2(i);
   end
   
   %figure
   subplot(2,1,2)
   %plot(temphist2)
   plot(densehist, 'LineWidth', 1.5, 'Color', 'Black')
   ylim([0 10]);
   xlabel('Time (sec)');
   ylabel('Number of spikes');
   title('ROI in Radius 130 \mum Spikes');
   
   
   
