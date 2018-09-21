clear;clf;
sc=1.8e-5; %scaling factor to offset the curves
%% load the data
file = 'SUP_20121109.mat';%'COM_20121109.mat'; 
load(file);

%% plot the computational basis and the modulation
figure(1);
subplot(1,3,1);
plot(t,real(y_mod)+0.22*(0:3)'*ones(1,size(y_mod,2)));
xlabel('Time');
ylabel('Field amplitude')
subplot(1,3,2);
plot(f,P+sc*(0:3)'*ones(1,size(P,2)));
xlabel('Frequency');
ylabel('Power /a.u.')
subplot(1,3,3);
plot(f,Pd+sc*(0:3)'*ones(1,size(Pd,2)));
xlabel('Frequency');
ylabel('Power /a.u.')

%% plot the field (real+imag) of the modulation
figure(2);
%plot(f,Py2+4e-1*(3:-1:0)'*ones(1,size(Py2,2)));
subplot(1,2,1);
plot(f,real(Y2+sqrt(sc)*(0:3)'*ones(1,size(Y2,2))));
title('Real part');
xlabel('Frequency');
ylabel('Field/a.u.');
subplot(1,2,2);
plot(f,imag(Y2+i*sqrt(sc)*(0:3)'*ones(1,size(Y2,2))));
title('Imaginary part');
xlabel('Frequency');
ylabel('Field/a.u.');

%% plot the modulation in the time domain
figure(3);
plot(t,abs(y_gs)+0.22*(0:3)'*ones(1,size(y_gs,2)));
%% plot the modulation basis and error matrix
figure(4);
subplot(1,2,1);
bar3(vec);
if strcmp(file(1:3),'SUP')
    title('Superposition basis');
else
    title('Compuational basis');
end
xlabel('Basis indices');
subplot(1,2,2);
%bar3(Pyx);
sanePColor(Pyx);colorbar;
title('Error matrix');
xlabel('Modulation indices');
ylabel('Demodulation indices');
