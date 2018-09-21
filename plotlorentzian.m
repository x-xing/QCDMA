N = 5e2;                            %number of steps
t_limit = 400;                      %time axis limit
t=linspace(-t_limit, t_limit, N);   %time axis in ns.
tau_p = 50;                         %photon lifetime.
param = [1, 1, 20];                 %[center, width] of a Lorentzian distribution for the filter.
n_fourier = 2^12;                   %dimension of fourier transform.

%Lorentzian
y = exp(-abs(t)/tau_p)/2/tau_p; % field integral is normalized.
%Gaussian
% y = exp(-t.^2/tau_p^2);
figure(1);
subplot(1,2,1);plot(t,y);title('Single photon field in time domain');xlabel('time (ns)')
Ys = fftshift(fft(y,n_fourier));
Pyy = Ys.* conj(Ys) /n_fourier;
f = ((-n_fourier/2):(n_fourier/2-1))/n_fourier;
zoomin = 10;
indzoom =round(((1-1/zoomin)/2*length(f)):((1+1/zoomin)/2*length(f)));
x1 = f(indzoom);
y1 = Pyy(indzoom);
subplot(1,2,2);plot(x1,[y1;fptransm(x1,param).*y1]');
title('Power spectrum of the single photon pulse');xlabel('frequency (GHz)');legend('Photon','Filter');
