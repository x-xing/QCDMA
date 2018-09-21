clear;
load('ErsLoss_20121109.mat');
dim = length(Py0x);
err = sum(Py0x(1:dim-1))/(Py0x(end)+sum(Py0x(1:dim-1)));
plot(err);