function Ixy = cha_cap_asym(param,sup_par)
%calculates the channel capacity. 
%param = [tau_EOM n_mod filt_width]
%sup_par = [sup_flag mod_depth(a)], par for superposition.

%param = [1 16 12]; sup_par=[1 0];

% Parameter initialization
tau_m = param(1);                          %modulation speed. minimum time to applied Vpi.
n_mod = param(2);
filt_w = param(3);
sup_flag = sup_par(1);
a = sup_par(2);                     %Amplitude modulation with depth a.
b = sqrt(1-a^2);
ab = [a; b];

dnu = 10e-3;                        %10MHz bandwidth of the photon
par_ph = [0, 2*pi*dnu/2];                 %[center, width] for the lorentzian of the photon.
f = (-1.2/tau_m):(dnu/20):(1.2/tau_m); 
af = f*2*pi;    %angular frequency, omega.

%par_ft = [peak transmission, FSR (GHz), finesse] for the filter.
par_ft = [1, 1.50, 1.5*2*pi/(filt_w*par_ph(2)*2)]; 
% ph_f = Ldistr(af,par_ph);    %the photon spectrum, not normalized
% filt_f = fptransm(f,par_ft);        %intensity transfer function
% figure(1);plot(f,ph_f,f,filt_f);

%eom speed determines the frequency shift
fshift = linspace(0,1/tau_m,n_mod)';
afshift = 2*pi*fshift;
onesSh = ones(size(afshift));
onesF  = ones(size(af));
Pyx = zeros(n_mod,n_mod);

if sup_flag
    %%%%%%%% superposistions %%%%%%%%%%
    for sp_v1 = 2%1:n_mod                 %basis to superimpose
        %mod_sp = b*Ldistr(onesSh*af-afshift*onesF,par_ph)...
        %    + a*Ldistr(onesSh*af-(afshift/2+onesSh*afshift(sp_v1)/2)*onesF,par_ph);
            
        Omg = onesSh*afshift(sp_v1);
        t_lim = 10/dnu;
        t = 0:tau_m/2.4:t_lim; %the sampling rate is 2.4/tau_m, or 2.2GHz here. Nyquist, 2/tau_m.
        Nt = numel(t);
        Nfft = 2^nextpow2(2*Nt);
        %note the 2pi factor, since omg = 2pi*f.
        %superposition of the type: a1+b2, b1-a2, a3+b4, b3-a4, etc.
        indodd = 1:2:n_mod;
        indeve = 2:2:n_mod;
        Omg1 = afshift(reshape(ones(2,1)*indodd,[],1));
        Omg2 = afshift(reshape(ones(2,1)*indeve,[],1));
        alt1 = (-1).^((0:n_mod-1)');
        coef1 = ab((alt1+3)/2);
        coef2 = ab((3-alt1)/2).*alt1;
        mod_sp_t = onesSh*exp(-abs(t-t_lim/2)*par_ph(2)).*(coef1*ones(1,Nt).*exp(i*Omg1*t)+coef2*ones(1,Nt).*exp(i*Omg2*t));
        %t-domain of the amplitude modulation only. check for resolution.
        % figure(3);plot(t,abs(mod_sp_t)+(1:n_mod)'*ones(size(t)));
        dm_t = (b+a*exp(-i*Omg*t))/sqrt(1+2*a*b);
        mod_spdm_t = mod_sp_t.*dm_t;
        mod_spdm =fftshift(fft(mod_spdm_t,Nfft,2),2);

        Fw = Nt/t_lim; %Frequency window, 1/t_step
        f1=linspace(-Fw/2,Fw/2,Nfft);
        onesF1 = ones(size(f1));
        mod_spdm = mod_spdm./(sqrt(sum(abs(mod_spdm).^2,2))*onesF1);  %normalization.
        mod_spdm = mod_spdm/sqrt(1+2*a*b); %losses.
        %time-domain(field) and power spectrum of the demodulation
%         figure(4);
%         subplot(1,2,1);plot(t,abs(mod_spdm_t)+(1:n_mod)'*ones(1,Nt));title('|amplitude|');
%         subplot(1,2,2);plot(f1,abs(mod_spdm).^2+0.01*(0:n_mod-1)'*onesF1);title('spectrum');
        det_sp = abs(mod_spdm).^2.*(onesSh*fptransm(f1-min(fshift),par_ft));
        Py0x = sum(det_sp,2);
    end
else
    mod = Ldistr(onesSh*af-afshift*onesF,par_ph);
    mod = mod./(sqrt(sum(abs(mod).^2,2)*onesF));   %normalization.
    %figure(2);plot(f,abs(mod).^2+0.01*(0:n_mod-1)'*onesF);
    %%power after the filter.
    det = abs(mod).^2.*(onesSh*fptransm(f-min(fshift),par_ft));
    %inter symbol interference. Probility of y0 given x(x0,x1,...,xn).
    Py0x = sum(det,2);
end

%%%% Note, by varying sp_v1, the Py0x value remains the same
%%%% thus, I used simiar shortcut to get the full Py0x.
%Probility of y(y0,y1,...,yn) given x(x0,x1,...,xn).
v0 = [Py0x;Py0x(n_mod:-1:2)];
for ind = 1:n_mod
    temp = circshift(v0,ind-1);
    Pyx(:,ind) = temp(1:n_mod);
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

