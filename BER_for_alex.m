load('BER_20120606.mat');
figure(1);
plot(err_lp');
legend('100','50','25','10');
xlabel('Hilbert space dimension');ylabel('BER');
figure(2);
plot(err_pf');
legend('100','50','25','10');
