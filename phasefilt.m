function phase1 = phasefilt(phase_pro, t,tau_m)
%model the finite speed of eom.

filt_n = round(tau_m*length(t)/max(t));               %average of the phase, not the field
filt_n = floor(filt_n/2)*2+1; %odd number of filter points, to be symmetric.
filt = 1/(sqrt(2*pi)*filt_n/4)*exp(-((1:filt_n)-round(filt_n/2)).^2/(2*(filt_n/4)^2));
filt=filt./sum(filt);
phase1 = [phase_pro(:,1)*ones(1,filt_n) phase_pro];  %padding the first filt_n datapoints of phase_pro.
phase1=filter(filt,1,phase1,[],2);         %moving average filtering.
phase1(:,1:filt_n)=[];                        %remove the padding.
