%%
%
%trains{1} = xspikesMEA(:,[1]);
%for i = 1:size(xspikesCalciumCells, 1)
%    trains{i+1} = xspikesCalciumCells{i}(:,[1]);
%end

%graphTrains(trains, 300);

%%
%
function graphTrains(TRAINS, duration, labels, offset)
    RASTERFIG = figure('WindowState','maximized');
    hold on
    %Graph top-down
    row = size(TRAINS,2);
    for i = 1:size(TRAINS,2)
        for j = 1:size(TRAINS{i},1)
            if size(TRAINS{i},2) > 0
                line([TRAINS{i}(j), TRAINS{i}(j)], [row , row+1], 'LineWidth', 2, 'Color', 'Black'); 
            end
        end
        %text(-offset,row+0.1,labels(i));
        text(-offset,row+.5,labels(i));
        row = row - 1;
    end
    set(gca,'TickDir','out') % draw the tick marks on the outside
    set(gca,'YTick', []) % don't draw y-axis ticks
    set(gca,'Color',get(gcf,'Color')) % match figure background
    %set(gca,'YColor',get(gcf,'Color')) % hide the y axis
    %set(gca,'XColor',get(gcf,'Color')) % hide the x axis
end