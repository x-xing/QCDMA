clear;clf;
load('ERSLoss20121113.mat');

semilogy(los,errloss');
%save('ErsLoss_20121109','Py0x');
legend('d=2','d=6','d=10','d=14','d=18')
xlabel('Loss')
ylabel('ERS')