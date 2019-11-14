function [sum_Spike_Count,mean_Spike_Rate,std_Spike_Count] = Spike_Counter(Spikes)

Spike_Count = sum(Spikes);
sum_Spike_Count = sum(Spike_Count);

mean_Spike_Count = mean(Spike_Count);
std_Spike_Count = std(Spike_Count)/size(Spikes,2);


Spike_Rate = (Spike_Count/size(Spikes,2));

mean_Spike_Rate = mean(Spike_Rate);