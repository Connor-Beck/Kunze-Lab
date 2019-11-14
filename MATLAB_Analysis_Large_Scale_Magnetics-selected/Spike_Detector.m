function Spikes = Spike_Detector(dF,std_threshold,static_threshold)

%dF function is the first derivative where
%dF(ROI_index,Time,Multiple_Videos)

Spikes = zeros(size(dF,1),size(dF,2));

dev = std_threshold*std(std(dF));
avghx = mean(mean(dF));
for j = 1:size(dF,1)
    for k = 1:size(dF,2)
        droploc=1;
        if k>droploc && sum(Spikes(j,1:k-1))==0
            if dF(j,k) >= dev + avghx && dF(j,k) >= static_threshold
                Spikes(j,k) = 1;
                location=k; %point at which dF surpasses the threshold
                for l=location:size(dF,2)
                    if dF(j,l) < 0
                        droploc = l; %time at which the cell decreases in intensity
                        break
                    end
                end
            end
        end
    end
end
