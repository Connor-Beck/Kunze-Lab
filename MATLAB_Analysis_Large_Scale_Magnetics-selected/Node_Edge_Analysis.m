
function [NumNodes,NumEdges] = Node_Edge_Analysis(ROIbaseH0,Connected_ROI)

Nodes = zeros(size(ROIbaseH0,1),length(Connected_ROI));
NumEdges = zeros(length(Connected_ROI),1);
NumNodes = zeros(length(Connected_ROI),1);
for i = 1:length(Connected_ROI)
    ROIlist = Connected_ROI{i};
    NumEdges(i) = length(ROIlist);
    for j = 1:length(ROIlist)
        C1 = ROIlist(j,1);
        C2 = ROIlist(j,2);
        if Nodes(C1,i) == 0
            Nodes(C1,i) = 1;  
        end
        if Nodes(C2,i) == 0
            Nodes(C2,i) = 1;
        end
    end
    NumNodes(i) = sum(Nodes(:,i));
end

subplot(1,2,1), plot(linspace(1,length(Connected_ROI),length(Connected_ROI)),NumNodes), ylim([0,1.1*max(NumNodes)])
subplot(1,2,2), plot(linspace(1,length(Connected_ROI),length(Connected_ROI)),NumEdges), ylim([0,1.1*max(NumEdges)])