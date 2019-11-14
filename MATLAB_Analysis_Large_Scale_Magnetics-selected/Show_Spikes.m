function spikeplotter = Show_Spikes(Spikes)

[h,w,d] = size(Spikes);
spikeplotter = zeros(h,w,d);
for i = 1:d
    for j = 1:h
        for k = 1:w
            if Spikes(j,k,i) == 1
                spikeplotter(j,k,i) = Spikes(j,k,i)+j;
            else
                spikeplotter(j,k,i) = Spikes(j,k,i)-10;
            end
        end
    end
end




for i = 1:d
    for j = 1:h
        subplot(1,d,i), scatter(linspace(0,w,w),spikeplotter(j,:,i),'k','.')
        hold on
        ylim ([0,h])
    end
end
