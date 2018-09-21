%calculates the BER for linear phase flip.

% Parameter initialization
filt_w = 1;
dnu = 10e-3;                        %10MHz bandwidth of the photon
tau_p = 1/dnu;                 
dim = 64;
err_pf=zeros(4,dim-1);
tm = [10,5,2,1];
t_lim = 20/dnu;

for ind = 1:4
    tau_m = tm(ind);
    t = 0:tau_m/2.4:t_lim; %the sampling rate is 2.4/tau_m, or 2.2GHz here. Nyquist, 2/tau_m.
    Nt = numel(t);
    Nfft = 2^nextpow2(2*Nt);
    y = exp(-abs(t-t_lim/2)*dnu*5)/2/tau_p;
    phase_pro = phapro(dim,t-t_lim/2,tau_m,tau_p,'pf');
    phase_pro = phasefilt(phase_pro, t,tau_m);
    y_mod = (ones(dim,1)*y).*exp(i*phase_pro);
    ys = fftshift(fft(y_mod,Nfft,2),2);
    ys = ys./(sqrt(sum(abs(ys).^2,2)*ones(1,Nfft)));   %normalization.
    f = linspace(-1,1,Nfft)/tau_m*1.2;
    Py = abs(ys).^2;
    %par_ft = [peak transmission, FSR (GHz), finesse] for the filter.
    par_ft = [1, 1.50, 1.5/(filt_w*dnu)];
    det = Py.*(ones(dim,1)*fptransm(f,par_ft));
    Py0x = sum(det,2);
    Py0x = sort(Py0x,'ascend');
    for n_mod = 2:dim
        err_pf(ind,n_mod-1) = mean(Py0x(1:n_mod-1))/(Py0x(end)+mean(Py0x(1:n_mod-1)));
    end
    disp([ind 4]);
end
%plot(f,Py(1:4:16,:));
plot(err_pf');

