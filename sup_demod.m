a=0.8;
N=2.4e4;2^nextpow2(1e4);
t=linspace(0,2e4,N);
Omg = 0.0741;

dm = (1-a)./(1+a*cos(Omg*t));%.*exp(-i*Omg*t);

subplot(2,1,1);plot(t,dm);
fdm = fftshift(fft(dm,10*N));
Pfdm = fdm.* conj(fdm) /N;
f = N/100*(1:10*N)/N; 
subplot(2,1,2);plot(f,Pfdm);

