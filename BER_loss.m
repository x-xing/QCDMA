
%calculates the BER for linear phase ramp.

% Parameter initialization
filt_w = 1;
darkratio = 1e-5;
dnu = 10e-3;                        %10MHz bandwidth of the photon
dim = 20;
modlist = 2:4:20;
errloss=zeros(10,length(modlist));
tm = [1, 2, 5, 10];
los = 1-logspace(-3,0,20);
for ind = 1:length(los)
    losind = los(ind);
    tau_m = tm(1);
    par_ph = [0, 2*pi*dnu/2];                 %[center, width] for the lorentzian of the photon.
    f = (-1.2/tau_m):(dnu/20):(1.2/tau_m);
    af = f*2*pi;    %angular frequency, omega.
    %par_ft = [peak transmission, FSR (GHz), finesse] for the filter.
    par_ft = [1, 15, 15*2*pi/(filt_w*par_ph(2)*2)];
    % ph_f = Ldistr(af,par_ph);    %the photon spectrum, not normalized
    % filt_f = fptransm(f,par_ft);        %intensity transfer function
    % figure(1);plot(f,ph_f,f,filt_f);
    for ind2=1:length(modlist)
        n_mod=modlist(ind2);
        dark = darkratio*n_mod;
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
        Py0x = Py0x*(1-losind) + dark;
        errloss(ind,ind2) = sum(Py0x(1:n_mod-1))/(Py0x(end)+sum(Py0x(1:n_mod-1)));
        %was the following: incorrect.
        %err_lp(ind,n_mod-1) = mean(Py0x(1:n_mod-1))/(Py0x(end)+mean(Py0x(1:n_mod-1)));
    end
end
%%
semilogy(los,errloss');
%save('ErsLoss_20121109','Py0x');
legend('d=2','d=6','d=10','d=14','d=18')
xlabel('Loss')
ylabel('ERS')