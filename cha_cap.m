%Note: obselete, fft parameter not optimized. See cha_cap_asym.m
function Ixy = cha_cap(param,sup_par)
%calculates the channel capacity. 
%param = [tau_EOM n_mod filt_width]
%sup_par = [sup_flag mod_depth(a)], par for superposition.

%param = [1 8 4]; sup_par=[1 0.1];
% Parameter initialization
tau_m = param(1);                          %modulation speed. minimum time to applied Vpi.
n_mod = param(2);
filt_w = param(3);
sup_flag = sup_par(1);
a = sup_par(2);                     %Amplitude modulation with depth a.


tau_p = 50;                         %photon lifetime.
f = linspace(-1.2/tau_m,1.2/tau_m,24e3/tau_m); % 0.1MHz resolution
par_ph = [0, 5e-3];                 %[center, width] for the lorentzian of the photon.
%par_ft = [peak transmission, FSR (GHz), finesse] for the filter.
par_ft = [1, 1.50, 1.5/(filt_w*par_ph(2)*2)]; 
% ph_f = Ldistr(f,par_ph);    %the photon spectrum, not normalized
% filt_f = fptransm(f,par_ft);        %intensity transfer function
% figure(1);plot(f,ph_f,f,filt_f);

%eom speed determines the frequency shift
fshift = linspace(0,1/tau_m,n_mod)';
onesSh = ones(size(fshift));
onesF  = ones(size(f));
Pyx = zeros(n_mod,n_mod);

if sup_flag
    %%%%%%%% superposistions %%%%%%%%%%
    for sp_v1 = 1:n_mod                 %basis to superimpose
        mod_sp = a/2*Ldistr(onesSh*f-fshift*onesF,par_ph)...
            + Ldistr(onesSh*f-(fshift/2+onesSh*fshift(sp_v1)/2)*onesF,par_ph)...
            + a/2*onesSh*Ldistr(f-fshift(sp_v1)*onesF,par_ph);

        Omg = fshift/2+onesSh*fshift(sp_v1)/2;
        t_lim = 20*tau_p;
        %make sure to have at least 10 points in the smallest cycle.
        N = 2^nextpow2(40*t_lim/tau_m);
        t = linspace(-t_lim,t_lim,N);
        %note the 2pi factor, since omg = 2pi*f.
        mod_sp_t = onesSh*exp(-abs(t)*par_ph(2)*2*pi).*exp(i*Omg*t).*(1+a*cos(Omg*t));
        %t-domain of the amplitude modulation only. check for resolution.
        %figure(3);plot(t,abs(mod_sp_t)+(1:n_mod)'*ones(size(t)));
        Omgd = Omg(sp_v1)*onesSh;
        dm_t = (1-a)./(1+a*cos(Omgd*t)).*exp(-i*Omgd*t);
        mod_spdm_t = mod_sp_t.*dm_t;
        mod_spdm =fftshift(fft(mod_spdm_t,N,2),2);
        %f-domain axis scaling. 0.1MHz.
        f1 = 1e5/N/(t(2)-t(1))/2*linspace(-1,1,N);
        onesF1 = ones(size(f1));
        mod_spdm = mod_spdm./(sqrt(sum(abs(mod_spdm).^2,2))*onesF1);  %normalization.
        mod_spdm = mod_spdm*(1-a)/sqrt(1+a^2/2); %losses.
        %time-domain(field) and power spectrum of the demodulation
        % figure(4);
        % subplot(1,2,1);plot(t,abs(mod_spdm_t)+(1:n_mod)'*onesT);
        % subplot(1,2,2);plot(f1,abs(mod_spdm).^2+0.25*(1:n_mod)'*onesF1);
        det_sp = abs(mod_spdm).^2.*(onesSh*fptransm(f1-min(fshift),par_ft));
        Pyx(:,sp_v1) = sum(det_sp,2);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    mod = Ldistr(onesSh*f-fshift*onesF,par_ph);
    mod = mod./(sqrt(sum(abs(mod).^2,2)*onesF));   %normalization.
    %figure(2);plot(f,abs(mod).^2+0.01*(1:n_mod)'*onesF);
    %%power after the filter.
    det = abs(mod).^2.*(onesSh*fptransm(f-min(fshift),par_ft));
    %inter symbol interference. Probility of y0 given x(x0,x1,...,xn).
    Py0x = sum(det,2);
    %Probility of y(y0,y1,...,yn) given x(x0,x1,...,xn).
    v0 = [Py0x;Py0x(n_mod:-1:2)];
    for ind = 1:n_mod
        temp = circshift(v0,ind-1);
        Pyx(:,ind) = temp(1:n_mod);
    end
end
%mutual information in this basis, assuming P(xi)=1/dim
%considering losses, P(yi)<1/dim = sum_j(P(yi|xj);
%Pyi: take into consideration the loss and noise.
%Losses: used for norm, but no entropy. 
%Noise:used for both norm and entropy
Pnorm = sum(Pyx)+1-max(Pyx);
Pyi = sum(Pyx)./Pnorm/n_mod;
Pxi = 1/n_mod*ones(size(Pyi));
Pyx = Pyx./(ones(length(Pyx),1)*Pnorm);
Ixy = mutinfo_c(Pxi,Pyi,Pyx);

