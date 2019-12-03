%%
%Get the EOI plot for a specific timeframe
startSec = 1;
endSec = 300;
EOItest = EOI(startSec*25000:endSec*25000);
clear startSec
clear endSec
plot(EOItest)
hline1 = refline(0,-4.5*EOIstd);
hline2 = refline(0,-5*EOIstd);
hline3 = refline(0,-5.5*EOIstd);
hline1.Color = 'r';
hline2.Color = 'r';
hline3.Color = 'r';

%%
