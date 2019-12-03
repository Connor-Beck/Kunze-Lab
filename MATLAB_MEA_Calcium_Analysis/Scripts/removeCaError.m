%% TEST
clear trains
for i = 1:size(xspikesCalciumCells, 1)
    trains{i} = xspikesCalciumCells{i}(:,[1]);
end
TEST_CLEANED = removeCaError(trains, 0.10);
clear labels
for i = 1:size(xspikesCalciumCells, 1)
    labels(i) = "";%horzcat("ROI " + num2str(i)); 
end
graphTrains(TEST_CLEANED, duration, labels, 5);
%%
test = TEST_CLEANED';
xspikesCalciumCells = test;
%for i = 1:size(xspikesCalciumCells, 1)
%    xspikesCalciumCells{i}(:,[1]) = TEST_CLEANED{i}
%end

%% RemoveCaError
%Removes noise for Ca spike trains
%Currently targets:
%           -Detection of most regions firing at the same time

function [CLEAN_TRAINS] = removeCaError(TRAINS_IN, min_percent_freq)
%Make an array for each specific Ca Time of fire. If new, place it
%inside. If duplicate, increment its count. Very inefficient, but
%should be okay. If the percentage of ROI firing at that timestamp
%exceeds the parameter, remove all of them.
%min_percent_freq is > 0 and <= 1
    CLEAN_TRAINS = TRAINS_IN;
    %COUNT THE FREQUENCIES
    num_keys = 0; %The number of fire timestamps in the table
    for i = 1:size(CLEAN_TRAINS,2) %each train
        for j = 1:size(CLEAN_TRAINS{i},1) %each spike
            found = 0;
            for k = 1:num_keys %search the table
                if timestamps(k, 1) == CLEAN_TRAINS{i}(j)
                    timestamps(k, 2) = timestamps(k, 2) + 1;
                    found = 1;
                end
            end
            if found == 0 %not found so add to table
                num_keys = num_keys + 1;
                timestamps(num_keys, 1) = CLEAN_TRAINS{i}(j);
                timestamps(num_keys, 2) = 1;
            end                
        end
    end
    
    %DELETE NECESSARY FIRES
    cutoff_freq = min_percent_freq * size(CLEAN_TRAINS,2);
    for k = 1:num_keys %iterate through the table
        if timestamps(k, 2) > cutoff_freq %if the fire time surpasses the cutoff, then find all occurences
            for i = 1:size(CLEAN_TRAINS,2) %each train
                j = 1;
                while (j <= size(CLEAN_TRAINS{i},1)) && (size(CLEAN_TRAINS{i},2) > 0)    %each spike
                    if CLEAN_TRAINS{i}(j) == timestamps(k, 1)
                       CLEAN_TRAINS{i}(j) = [];
                       j = j - 1;
                    end
                    j = j + 1; 
                end
            end
        end
    end
end