
%calculates the BER for linear phase ramp.

% Parameter initialization
filt_w = 1;
dnu = 10e-3;                        %10MHz bandwidth of the photon
dim = 64;
err_lp=zeros(4,dim-1);
tm = [1, 2, 5, 10];
for ind = 1
    tau_m = tm(ind);
    par_ph = [0, 2*pi*dnu/2];                 %[center, width] for the lorentzian of the photon.
    f = (-1.2/tau_m):(dnu/20):(1.2/tau_m);
    af = f*2*pi;    %angular frequency, omega.
    %par_ft = [peak transmission, FSR (GHz), finesse] for the filter.
    par_ft = [1, 15, 15*2*pi/(filt_w*par_ph(2)*2)];
    % ph_f = Ldistr(af,par_ph);    %the photon spectrum, not normalized
    % filt_f = fptransm(f,par_ft);        %intensity transfer function
    % figure(1);plot(f,ph_f,f,filt_f);
    for n_mod = 2:dim
        %eom speed determines the frequency shift
        fshift = linspace(0,1/tau_m,n_mod)';
        afshift = 2*pi*fshift;
        onesSh = ones(size(afshift));
        onesF  = ones(size(af));
        mod = Ldistr(onesSh*af-afshift*onesF,par_ph);
        mod = mod./(sqrt(sum(abs(mod).^2,2)*onesF));   %normalization.
        %%power after the filter.
        det = abs(mod).^2.*(onesSh*fptransm(f-min(fshift),par_ft));
        %inter symbol interference. Probility of y0 given x(x0,x1,...,xn).
        Py0x = sum(det,2);
        Py0x = Py0x(end:-1:1);
        disp([ind n_mod dim]);
        err_lp(ind,n_mod-1) = sum(Py0x(1:n_mod-1))/(Py0x(end)+sum(Py0x(1:n_mod-1)));
        %was the following: incorrect.
        %err_lp(ind,n_mod-1) = mean(Py0x(1:n_mod-1))/(Py0x(end)+mean(Py0x(1:n_mod-1)));
    end
end
plot(err_lp');
legend('100','50','20','10');
title('linear phase ramp');
%save('ErsLoss_20121109','Py0x');