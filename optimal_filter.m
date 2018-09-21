f= -100:0.1:100;
param1 = [1, 200, 100];            %[peak transmission, FSR (GHz), finesse] for the filter.
param2 = [1, 200, 100];            %lorentzian for the photon.
f0 = 50;

tp = zeros(100,2);
for ind = 1:length(tp)
    param1(3) = param1(2)/ind*25;
    eff = fptransm(f,param1).*fptransm(f,param2);
    err = fptransm(f-f0,param1).*fptransm(f,param2);
    tp(ind,:) = [sum(eff) sum(err)];
end


subplot(3,1,1);plot(tp(:,1));title('eff');
subplot(3,1,2);plot(tp(:,2));title('err');
subplot(3,1,3);plot(tp(:,1).^2./tp(:,2));title('eff^2/err');xlabel('filter width (25*photon)');