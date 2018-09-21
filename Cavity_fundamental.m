p1=pwd;
p1=[p1(1:strfind(p1, 'Dropbox')-1) 'Dropbox\Data\QCDMA Filter Cavity'];

f0ch12=load([p1 '\F0000CH12.CSV']);
f1ch12=load([p1 '\F0001CH12.CSV']);

subplot(2,2,1);
plot(f0ch12(:,1),f0ch12(:,2));
subplot(2,2,3);
plot(f0ch12(:,1),f0ch12(:,4));
subplot(2,2,2);
plot(f1ch12(:,1),f1ch12(:,2));
subplot(2,2,4);
plot(f1ch12(:,1),f1ch12(:,4));
