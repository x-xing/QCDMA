clear all;
load('mod_data.mat');
%variable: t, y_mod (time domain); f, Py( frequency domain)

subplot(1,2,1);
plot(t,real(y_mod(1:3:16,:))+1e-2*(0:5)'*ones(1,size(y_mod,2)));
subplot(1,2,2);
plot(f,Py(1:3:16,:)+2.5e-2*(0:5)'*ones(1,size(Py,2)));