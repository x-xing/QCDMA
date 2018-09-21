close all;
clear;
load('BER_20120925.mat');
figure(1);
hold on;
plot(err_pf_d16e1,'bs-');
plot(err_pf_d32e1,'gx-');
plot(err_pf_d64e1,'rd-');
hold off;
legend('n=16','n=32','n=64');
xlabel('Hilbert space dimension');ylabel('ERS');title('ERS of phase flip scheme: EOM/photon = 100')
figure(4);
hold on;
plot(err_pf_d16e5,'bs-');
plot(err_pf_d32e5,'gx-');
plot(err_pf_d64e5,'rd-');
hold off;
legend('n=16','n=32','n=64');
xlabel('Hilbert space dimension');ylabel('ERS');title('ERS of phase flip scheme: EOM/photon = 20')

figure(2);
pcol = ['bs-';'gx-';'rd-';'k^-'];
hold on;
for ind = 1:size(err_lp,1)
    plot(err_lp(ind,:),pcol(ind,:));
end
hold off;
legend('N=10','N=20','N=50','N=100');
title('ERS of linear phase ramp scheme');
xlabel('Hilbert space dimension');ylabel('ERS');
figure(3);
hdim = ones(4,1)*(1:63)./([10;20;50;100]*ones(1,63));
newerr=err_lp';
hold on;
for ind=1:4
    plot(log(hdim(ind,:)),log(newerr(:,ind)),pcol(ind,:));
end
hold off;
title('ERS of linear phase ramp scheme');
xlabel('log_e((d-1)/N)');
ylabel('log_e(ERS)');    
legend('N=10','N=20','N=50','N=100');    
