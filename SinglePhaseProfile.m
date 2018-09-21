a=(phase_pro(2,:));
af=filter(filt,1,[ones(1,filt_n)*a(1) a]);
plot(af)