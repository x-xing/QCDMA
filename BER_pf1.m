%calculates the BER for linear phase flip.

%% Parameter initialization
filt_w = 1;
dnu = 10e-3;                        %10MHz bandwidth of the photon
tau_p = 1/dnu;                 
dim = 32;
err_pf=zeros(4,dim-1);
tm = [10,2,1];
t_lim = 20/dnu;
err = [];
clear i;
step = 5;

%% compare dimension to eom speed
% smallest division -ln(1-1/dim) < tm/tp ??

%% Generating the phase profile and calculate Pyx
for ind = 1:1
    tau_m = tm(ind);
    t = 0:tau_m/step:t_lim; %the sampling rate is 2.4/tau_m, or 2.4GHz here. Nyquist, 2/tau_m.
    Nt = numel(t);
    Nfft = 2^nextpow2(Nt);
    y = exp(-abs(t-t_lim/2)*dnu)/2/tau_p;
    disp('Generating phase profile...');
    tic;
    phase1 = phapro1(dim,t-t_lim/2,tau_m,tau_p,'pf'); %phapro1: ideal phase profile
    %Model the finite EOM speed. phase1 is already the product!!! 
    phase1=phasefilt(phase1, t,tau_m); 
    y_mod = (ones(dim,1)*y).*exp(i*phase1);
    toc;
    tic;
    disp('FFT...');
    ys = fftshift(fft(y_mod,Nfft,2),2);
    toc;
    ys = ys./(sqrt(sum(abs(ys).^2,2)*ones(1,Nfft)));   %normalization.
    f = linspace(-1,1,Nfft)/tau_m*step/2;
    Py = abs(ys).^2;
    %par_ft = [peak transmission, FSR (GHz), finesse] for the filter.
    par_ft = [1, 1.50, 1.5/(filt_w*dnu)];
    det = Py.*(ones(dim,1)*fptransm(f,par_ft));
    Py0x = sum(det,2); % Spectrum of the product of two profiles.
end
%generating the full matrix Pyx from product Py0x.
R = zeros(dim);
for i1 = 0:dim-1
    R(i1+1,:) = bitxor(0:dim-1, i1*ones(1,dim));
end
Pyx = Py0x(R+1);
% add the effect of EOM here???

figure(1);
pcolor(Pyx);colormap('hot');colorbar;
xlabel('Input basis state index');
ylabel('Output basis state index');
title(['Error matrix for phase flip encoding, N=' num2str(100/tm(1))]);
s1 = sprintf('EOM speed: %d',100/tm(ind));



%% Optimization by exhaustive search
tic;
if dim <= 16
    disp('Optimizing...');
    %[err indcell] = beropt(Pyx);
else
    disp('Skipping brute force calculation...');
end
toc;


%% Approximate optimization
figure(2);
err1 = bersimpopt(Pyx,tau_m/tau_p);
if err 
    plot([err;err1]');
else
    hold on; plot(err1);
end
