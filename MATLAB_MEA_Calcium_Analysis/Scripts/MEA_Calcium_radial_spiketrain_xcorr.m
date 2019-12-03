%% MEA Calcium Radial Spiketrain Cross Correlation
%@Andy Kirby
% Script designed to take in spike trains from MEA and Calcium data, create
% spiketrains based on radius from an electrode, and perform cross
% correlation on each with respect the the electrode spiketrain.
%
% Currently, pairs with the MEA Calcium Analysis script, as it directly
% uses variables from that script's workspace.

% Define radii to test
start_rad = 60;
max_rad = 800;
radii_interval = 5;

num_radii = (max_rad - start_rad) / radii_interval;
radii_to_test = zeros(num_radii, 1);%[60,71,90,150,350];
for i = 1:num_radii
    radii_to_test(i) = start_rad + (radii_interval*(i-1));
end

%Create array to store chi values
% chi_vals = zeros(size(radii_to_test, 1), 2);
% hist_bin_size = 15;
hist_bin_size = 1;

% Create an array to store cross correlation values
% correlations = zeros(length)

    
%bin_data holds the binned spikes matrix, possibly include import ability
%input or import ROI distances/radii
ROI_distances = getDistances(electrode_coord, ROIboundary, 0.426727);%[62,54,70,54,58,62,70,87,74,82,111,120,136,144,185,243,268,338,309];

%for radius in radius-to-test list
%Read radius from inputted radius-to-test list
for r = 1:size(radii_to_test,1) 
    
    cur_radius = radii_to_test(r); %in um
    %Create the matrix that each test will be run on
    hist_compare = zeros(size(bin_data,1)/hist_bin_size,3);
    %Load hist_compare with appropriate distributions
    for i = 1:size(hist_compare,1)
        for k = hist_bin_size*(i-1)+1:hist_bin_size*(i-1)+hist_bin_size
            %hist_compare(i,1) = hist_compare(i,1) + bin_data(k,1); %copy time col
            hist_compare(i,1) = bin_data(k,1); %copy time col
            hist_compare(i,2) = hist_compare(i,2) + bin_data(k,2); %copy MEA col
            %Add Ca ROI spikes only if they fall into radius
            for j = 3:size(bin_data,2)
                if ROI_distances(j-2) <= cur_radius
                    hist_compare(i,3) = hist_compare(i,3) + bin_data(k,j);
                end
            end
        end
    end
    
    % Apply Cross Correlation
    MEA_spiketrain = hist_compare(:,2);
    RADIUS_spiketrain = hist_compare(:,3);
    c = xcorr(MEA_spiketrain, RADIUS_spiketrain,0,'normalized');
    c = max(c);
    correlations(r,1) = radii_to_test(r);
    correlations(r,2) = c;
end

%% Graph the xcorr vs radius
figure;
plot(correlations(:,1), correlations(:,2));
title('MEA-Calcium Spiketrain Cross Correlation vs Radius');
xlabel('Radius (\mum)');
ylabel('XCORR Coefficient');

%%
%Write the data table to a file
%Start with the header
fid = fopen(horzcat('../Output/', 'correlation_vs_radius.csv'), 'w');
header={'RADIUS','XCORR'};
fprintf(fid,'%s,%s',header{:});
fprintf(fid,'\n');
fclose(fid);
%Write the data
dlmwrite(horzcat('../Output/', 'correlation_vs_radius.csv'), correlations, '-append');

