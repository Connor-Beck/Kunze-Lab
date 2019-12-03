function [xspikesMEA, tspikesMEA, EOI_copy, NCED] = ncedDetect(EOI, xspikesMEA, tspikesMEA, Fs)

    % ~~~~Edit:Andy Kirby~~~~
    % Utilize Normalized Cumulative Energy Difference (NCED) 
    % to detect spikes per Mtetwa and Smith (2006).
    % First, find the cumulative energies of the signal, Ecumulative, and the
    % total energy of the signal, Etot, by summing the squares of each 
    % recorded value in the signal.

    %Make a copy of the EOI that eliminates all positive values
    EOI_copy = EOI;
    for i = 1:size(EOI_copy)
        if EOI_copy(i) > 0
            EOI_copy(i) = 0;
        end
    end

    Ecumulative = zeros(size(EOI_copy,1),1);
    Etotal = 0;
    for i = 1:size(EOI_copy)
        Etotal = Etotal + EOI_copy(i)^2;
        Ecumulative(i) = Etotal;
    end
    % Normalize the cumulative energy
    Ecumulative = Ecumulative / Etotal;

    % Next, find the differences of the NCEs (NCED) that will be used to 
    % identify spikes.
    NCED = zeros(size(EOI_copy,1),1);
    for i = 2:size(EOI_copy)
        NCED(i) = Ecumulative(i) - Ecumulative(i-1);
    end
    % Account for dt, our timestep interval, which is related to the capture
    % rate, Fs.
    NCED = NCED * Fs;

    %spikeDETECTED = 0
    %for i = 1:size(EOI)
    %    if NCED(i) > 1
    %        spikeDETECTED = 1
    %    end
    %end

    %~~~~~~~~~~~~~~~~~~~~~~~~

    %%
    % EOIthresh is 1 due to NCED spikes occuring when the value is much larger
    % than 1
    EOIthresh = 1;
    % Initialize below threshold binary variable
    if NCED(1) < EOIthresh
        belowthresh = 1;
    else
        belowthresh = 0;
    end

    % Comb through EOI data. If value is below threshold and data was
    % previously above threshold, mark the time and correct the below threshold 
    % variable. If value is above threshold and data was previously below 
    % threshold, just correct the below threshold variable.
    for i = 1:size(EOI)
        if NCED(i) < EOIthresh && belowthresh == 0
            tspikesMEA(i) = 1;
            belowthresh = 1;
        elseif NCED(i) > EOIthresh && belowthresh == 1
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

end