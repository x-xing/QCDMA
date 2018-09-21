close all;
load('berf20110725.mat');

figure (1);

plot(ber_m0p1_type_2(1,:),ber_m0p1_type_2(2,:),ber_m0p2_type_2(1,:),ber_m0p2_type_2(2,:),...
ber_m0p3_type_2(1,:),ber_m0p3_type_2(2,:),ber_m0p4_type_2(1,:),ber_m0p4_type_2(2,:),ber_m0p5_type_2(1,:),ber_m0p5_type_2(2,:),...
ber_m0p6_type_2(1,:),ber_m0p6_type_2(2,:),ber_m0p7_type_2(1,:),ber_m0p7_type_2(2,:),ber_m0p8_type_2(1,:),ber_m0p8_type_2(2,:),...
ber_m0p9_type_2(1,:),ber_m0p9_type_2(2,:),ber_m1p0_type_2(1,:),ber_m1p0_type_2(2,:));

legend('0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1');
title('BER for phase flip for different EOM speed');
xlabel('dimension');
ylabel('BER');

figure (2);

plot(ber_m0p1_type_3(1,:),ber_m0p1_type_3(2,:),ber_m0p2_type_3(1,:),ber_m0p2_type_3(2,:),...
ber_m0p3_type_3(1,:),ber_m0p3_type_3(2,:),ber_m0p4_type_3(1,:),ber_m0p4_type_3(2,:),ber_m0p5_type_3(1,:),ber_m0p5_type_3(2,:),...
ber_m0p6_type_3(1,:),ber_m0p6_type_3(2,:),ber_m0p7_type_3(1,:),ber_m0p7_type_3(2,:),ber_m0p8_type_3(1,:),ber_m0p8_type_3(2,:),...
ber_m0p9_type_3(1,:),ber_m0p9_type_3(2,:),ber_m1p0_type_3(1,:),ber_m1p0_type_3(2,:));

legend('0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1');
title('BER for linear phase for different EOM speed');
xlabel('dimension');
ylabel('BER');

