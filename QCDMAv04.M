%modified from QCDMAv02. simplify the structure of calculation, simulation
%of actual modulation-demodulation process. 
% -- Feb. 08, 2012

%clear all;
% tpft = [];
close all;
%% Use qopath to add path.
p1=pwd;
cd([p1(1:strfind(p1, 'Dropbox')-1) 'Dropbox\MatlabCode\utility']);
qopath;
cd(p1);

%% Parameter initialization.
N = 8e3;                            %number of steps 4e4 for PF
t_limit = 200;                      %time axis limit
t=linspace(-t_limit, t_limit, N);   %time axis in ns.
tau_p = 100;                         %photon lifetime.
tau_m = 1;                         %modulation speed. minimum time to applied Vpi.
param = [1, 1.50, 1500];            %[peak transmission, FSR (GHz), finesse] for the filter.
plotsingle = 1;                     %parameter to plot a single dimension, or a range.
modtype = 'lp';                        %modulation type. 1:UCP, 2: PF, 3: LP
n_mod = 8;                          %dimension of the encoding. valid in single dimension cal.
n_mod_max =8;%round(1/(1-exp(-8*tau_m/tau_p))); %to make sure time bin is larger than eom speed.
n_mod_range = 2:n_mod_max;                 %range of dimension of the encoding. dimension is a variabhhle.
n_fourier = 2^nextpow2(N);                   %dimension of fourier transform.

%typical setting: as a function of dimension: 
%plotsingle = 0; 
%modtype = 'lp', (linear phase); n_mod_range = 2:15
%modtype = 'pf', (phase filp); n_mod_range = 2:50.


%% Plot the time and frequency domain of the single photon 
%Lorentzian, exponential in time.
y = exp(-abs(t)/tau_p)/2/tau_p; % field integral is normalized.
%Gaussian % y = exp(-t.^2/tau_p^2);
figure(1);
subplot(1,2,1);plot(t,y);title('Single photon field in time domain');xlabel('time (ns)')
Ys = fftshift(fft(y,n_fourier));
Pyy = Ys.* conj(Ys) /n_fourier;
f = N/(2*t_limit)*((-n_fourier/2):(n_fourier/2-1))/n_fourier; % see Matlab fft example for the pre-factor.
zoomin = 100;
indzoom =round(((1-1/zoomin)/2*length(f)):((1+1/zoomin)/2*length(f)));
x1 = f(indzoom);
y1 = Pyy(indzoom);
subplot(1,2,2);plot(x1,[y1;fptransm(x1,param)*max(y1)]');  %fptransm power transmission, here just plot to see.
title('Power spectrum of the single photon pulse');xlabel('frequency (GHz)');legend('Photon','Filter');

effavg = zeros(1,n_mod_max-1);
effmin = zeros(1,n_mod_max-1);
erravg = zeros(1,n_mod_max-1);
errmax = zeros(1,n_mod_max-1);



if plotsingle
    n_mod_range = n_mod;
end
n_mod_last = 2;

for n_mod = n_mod_range
    sprintf('%d of %d',n_mod,max(n_mod_range))

    if (n_mod < nextdim(n_mod)) && strcmp(modtype, 'pf') && ~plotsingle
        continue;
    else if strcmp(modtype, 'pf') && plotsingle
        n_mod = nextdim(n_mod);
        end
    end
    phase_pro = phapro(n_mod,t,tau_m,tau_p,modtype);
    y_mod = (ones(n_mod,1)*y).*exp(i*phase_pro);

    %%%%%%%%%%%
    %superpostion calculation, for linear ramp only. 
    %phase_pro = pi/tau_m/(dim-1)*(0:(dim-1))'*t;
    sp_flag = 0;
    sp_v1 = 2; %basis used for superpositions.
    Omg = pi/tau_m/(n_mod-1)*(sp_v1-1)/2;
    a = 0.5;
    %demodulation, with frequency shift.
    dm = (1-a)./(1+a*cos(Omg*t)).*exp(-i*Omg*t);
    dm = 1; %use this to see the modulated spectrum.
    mod_sp = exp(i*phase_pro)+exp(i*ones(n_mod,1)*phase_pro(sp_v1,:));
    %Amplitude modulation with depth a.
    mod_sp = (mod_sp*a/2 + exp(i*(.5*phase_pro+.5*ones(n_mod,1)*phase_pro(sp_v1,:))))/(1+a^2/2);
    y_mod = (ones(n_mod,1)*(y.*dm)).*mod_sp;
    %%%%%%%%%%%
    
    y_mod_eff = abs(y_mod);                      %this is for the efficiency calculation.

    Ys_mod = fftshift(fft(y_mod',n_fourier),1)';
    Ys_mod_eff = fftshift(fft(y_mod_eff',n_fourier),1)';
    Pyy_mod = Ys_mod.* conj(Ys_mod) / n_fourier;
    Pyy_mod_eff = Ys_mod_eff.* conj(Ys_mod_eff) / n_fourier;

    %Calculate the efficiency of detection with a cavity filter.
    err = sum(Pyy_mod.*(ones(n_mod,1)*fptransm(f,param)),2)/sum(Pyy);
    eff = sum(Pyy_mod_eff.*(ones(n_mod,1)*fptransm(f,param)),2)/sum(Pyy);
    %pick the minimum error rate encoding.
    [err1_st, ind_st]=sort(err,'ascend');
    
    err1 = err1_st(1:n_mod);
    eff1 = eff(ind_st(1:n_mod));
    y_mod = y_mod(ind_st(1:n_mod),:);
    Pyy_mod = Pyy_mod(ind_st(1:n_mod),:);
    %y_mod = y_mod(ind_st(1:4:n_mod),:);         %for the QCDMA paper
    %Pyy_mod = Pyy_mod(ind_st(1:4:n_mod),:);

    effavg(n_mod-1) = mean(eff1);
    effmin(n_mod-1) = min(eff1);
    erravg(n_mod-1) = mean(err1(1:(n_mod-1)));
    errmax(n_mod-1) = max(err1(1:(n_mod-1)));
    if strcmp(modtype, 'pf')
        for ind = n_mod_last:n_mod
            err1 = err1_st(1:ind);
            eff1 = eff(ind_st(1:ind));
            effavg(ind-1) = mean(eff1);
            effmin(ind-1) = min(eff1);
            erravg(ind-1) = mean(err1(1:(ind-1)));
            errmax(ind-1) = max(err1(1:(ind-1)));
            %ind = ind + 1;
        end
    end

    if plotsingle
        figure(2);
        subplot(1,2,1);plot(t,real(y_mod'+ones(length(t),1)*2e-2*(1:length(y_mod(:,1)))));
        title('modulated field (real part) in time domain');xlabel('time (ns)');
        zoom = 10;
        indzoom =round(((1-1/zoom)/2*length(f)):1:((1+1/zoom)/2*length(f)));
        subplot(1,2,2);plot(f(indzoom),Pyy_mod(:,indzoom)'+2.5e-2*ones(length(indzoom),1)*(0:(length(y_mod(:,1))-1)));
        hold on;
        plot(f(indzoom),max(Pyy_mod(1,indzoom))*fptransm(f(indzoom),param),'k-.');hold off;
        title('Power spectrum of the modulated pulse');xlabel('frequency (GHz)');

        figure(3);
        len1 = length(err1)-1;
        plot(err1(1:len1),'bo');
        set(gca,'XTickLabel',num2cell(ind_st(1:len1)));     
        xlabel('basis index');
        title([modtype ': eff_{avg}=' num2str(mean(eff1)*100)...
            '%, eff_{min}=' num2str(min(eff1)*100)...
            '%,   err_{avg}=' num2str(mean(err1(1:len1))*100)...
            '%, err_{max}=' num2str(max(err1(1:len1))*100) '%']);
        ylabel('error rate');
        
        %tpft = [tpft; mean(eff1)/(mean(eff1)+mean(err1(1:(length(err1)-1)))) param(2)/param(3)*1000];
    end
    n_mod_last = n_mod;
end

if ~plotsingle
    figure(2);
    subplot(1,2,1);plot(n_mod_range,effavg,'bo',n_mod_range,effmin,'rx');
    title('Average/Min efficiency as a function of dimension'); xlabel('Encoding Dimension');ylabel('Efficiency');
    legend('Average efficiency','Min efficiency');
    subplot(1,2,2);plot(n_mod_range,erravg,'bo',n_mod_range,errmax,'rx');
    title('Average/Max error rate as a function of dimension'); xlabel('Encoding Dimension');ylabel('error rate');
    legend('Average error rate','Max error rate');
    figure(3);
    plot(n_mod_range,erravg./(effavg+erravg),'bo');
    title('BER = Err_{avg}/all'); xlabel('Encoding Dimension');ylabel('Effective BER');

    vname = sprintf('ber_m%dp%d_type_%d',floor(tau_m),(tau_m-floor(tau_m))*10, modtype);
    v=genvarname(vname);
    eval([v  '=[n_mod_range; erravg./effavg];']);
    fname = sprintf('berf%s.mat',datestr(now, 'yyyymmdd'));

    f1=exist(fname);
    if f1
        save(fname,vname,'-append');
    else
        save(fname,vname);
    end
end

