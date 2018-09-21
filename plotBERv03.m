close all;
clear all;
load('berf20120102.mat');

figure (1);

plot(ber_m1p0_type_2(1,:),ber_m1p0_type_2(2,:),ber_m2p0_type_2(1,:),ber_m2p0_type_2(2,:),...
    ber_m5p0_type_2(1,:),ber_m5p0_type_2(2,:),ber_m10p0_type_2(1,:),ber_m10p0_type_2(2,:));

legend('50','25','10','5');
title('BER for phase flip for different EOM speed');
xlabel('dimension');
ylabel('BER');

figure (2);

plot(ber_m1p0_type_3(1,:),ber_m1p0_type_3(2,:),ber_m2p0_type_3(1,:),ber_m2p0_type_3(2,:),...
     ber_m5p0_type_3(1,:),ber_m5p0_type_3(2,:),ber_m10p0_type_3(1,:),ber_m10p0_type_3(2,:));

legend('50','25','10','5');
title('BER for linear phase for different EOM speed');
xlabel('dimension');
ylabel('BER');

