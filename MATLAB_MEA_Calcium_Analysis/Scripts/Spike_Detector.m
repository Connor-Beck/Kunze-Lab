function Spikes = Spike_Detector(DeltaFoverF)


Spikes = zeros(size(DeltaFoverF,1),size(DeltaFoverF,2),size(DeltaFoverF,3));

for i = 1:size(DeltaFoverF,3)
    dev = 5*std(std(DeltaFoverF(:,:,i)));
    for j = 1:size(DeltaFoverF,1)
        for k = 1:size(DeltaFoverF,2)
            if DeltaFoverF(j,k,i) >= dev
                Spikes(j,k,i) = 1;
            end
        end
    end
end
   
    