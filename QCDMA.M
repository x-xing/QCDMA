%outdated.

clear all;
close all;
%% Use qopath to add path.
p1=pwd;
cd([p1(1:strfind(p1, 'Dropbox')-1) 'Dropbox\MatlabCode\utility']);
qopath;
cd(p1);

%% Parameter initialization.
N = 500;                            %number of steps
t_limit = 400;                      %time axis limit
t=linspace(-t_limit, t_limit, N);   %time axis in ns.
tau_p = 50;                         %photon lifetime.
tau_m = 1;                          %modulation speed.
param = [0., 2e-3];                 %[center, width] of a Lorentzian distribution for the filter.
plotsingle = 1;                     %parameter to plot a single dimension, or a range.
modtype = 3;                        %modulation type. 1:UCP, 2: PF, 3: LP
n_mod = 2;                          %dimension of the encoding. valid in single dimension cal.
n_mod_range = 2:12;                 %range of dimension of the encoding. dimension is a variable.
n_fourier = 2^12;                   %dimension of fourier transform.

%typical setting: as a function of dimension: 
%plotsingle = 0; 
%modtype = 3, (linear phase); n_mod_range = 2:15
%modtype = 2, (phase filp); n_mod_range = 2:50.


%% Plot the time and frequency domain of the single photon 
%Lorentzian
y = exp(-abs(t)/tau_p)/2/tau_p; % field integral is normalized.
%Gaussian
% y = exp(-t.^2/tau_p^2);
figure(1);
subplot(1,2,1);plot(t,y);title('Single photon field in time domain');xlabel('time (ns)')
Ys = fftshift(fft(y,n_fourier));
Pyy = Ys.* conj(Ys) /n_fourier;
f = ((-n_fourier/2):(n_fourier/2-1))/n_fourier;
zoomin = 5;
indzoom =round(((1-1/zoomin)/2*length(f)):((1+1/zoomin)/2*length(f)));
x1 = f(indzoom);
y1 = Pyy(indzoom);
subplot(1,2,2);plot(x1,[y1;Ldistr(x1,param)*max(y1)]');
title('Power spectrum of the single photon pulse');xlabel('frequency (GHz)');legend('Photon','Filter');

effavg=[];
effmin=[];
erravg=[];
errmax=[];
if plotsingle
    n_mod_range = n_mod;
end

%% generate the phase profile
for n_mod =n_mod_range

    %%1. Unit-circle phase encoding (UCP)
    %%to-do
    if modtype ==1                 
        phase_mod = exp(i*2*pi*(1:n_mod)/n_mod);           %phase changes in steps.
    end
    
    %%2. Phase-flip encoding (PF)
    % Hadamard matrix (n x n) exists only when rem(n,4)=0. 
    % Matlab hadamard(n) function handles only the cases where n, n/12, or n/20 is a power of 2.
    if modtype ==2
        d1= min(find(n_mod./2.^(1:50)<=1));                 %the minimum number to satisfy 2.^d > n_mod.
        d2= min(find(n_mod/12./2.^(0:50)<=1));              %the minimum number to satisfy 12*2.^d > n_mod.
        d3= min(find(n_mod/20./2.^(0:50)<=1));              %the minimum number to satisfy 20*2.^d > n_mod.
        n_modax = min([2^d1; 12*2^(d2-1); 20*2^(d3-1)]); %increase the dimension as required by hadamard.
        hada_d=hadamard(n_mod);                               
        %phase_mod=hada_d(2,:);                             %alternating in every consective bins.
        phase_mod=hada_d;                                  %every possible phase modulation vector.
        %count the number of flips
        pdiff = diff(phase_mod,1,2);
        num_flip = sum(abs(pdiff)/2,2);
    end
    
    %%Calculate the phase change point for PF.
    if modtype == 2
        int_pul = int_exp(t,tau_p);                         %integrated the lorentzian from -inf to t.
        %%% find the indices of the phase change(flip) points. To verify, use int_pul(ind_t).
        ind_t = zeros(1,n_mod); 
        for ind_m = 1:n_modax
            [val,ind_t(ind_m)] = min(abs(int_pul-ind_m/n_mod)); 
        end
        phase_pro = zeros(n_mod, length(t));                %phase profile matrix.
        phase_pro(:, 1:ind_t(1)) = (phase_mod(:, ones(1,ind_t(1)))-1)*pi/2;
        for ind_m = 1:(n_mod-1)
            phase_pro(:, ind_t(ind_m):ind_t(ind_m+1)) = (phase_mod(:,(ind_m+1)*ones(1,ind_t(ind_m+1)-ind_t(ind_m)+1))-1)*pi/2;
        end
    end
    
    %%3. Linear phase encoding (LP)         
    %%to-do
    if modtype ==3
        N = 8e3;                                    %10 points in every pi encoding.
        t=linspace(-t_limit, t_limit, N);           %time axis in ns.
        y = exp(-abs(t)/tau_p)/2/tau_p;             %field integral is normalized.
        Ys = fftshift(fft(y,n_fourier));
        Pyy = Ys.* conj(Ys) / n_fourier;
        phase_pro = zeros(n_mod, length(t));         %phase profile matrix.
        phase_pro(1,:)=ones(1,length(t));
        phase_pro(2:n_mod,:)= pi/tau_m/n_mod.*(1:(n_mod-1))'*t; %%%%%%%%start from 0?????????????????????????????????????????
    end
    
    %% Model the finite EOM speed with a gaussian kernel.
    % Using the moving average method here, with a gaussian weights whose sig = tau_m.
    filt_n = round(tau_m*5*N/t_limit);                            % average of the phase, not the field
    filt = 1/(sqrt(2*pi)*tau_m)*exp(-((1:filt_n)-round(filt_n/2)).^2/(2*tau_m^2));
    phase_pro = [phase_pro(:,1)*ones(1,filt_n) phase_pro];        %padding the first filt_n datapoints of phase_pro.
    phase_pro=filter(filt,1,phase_pro,[],2);                      %moving average filtering.
    phase_pro(:,1:filt_n)=[];                                     %remove the padding.
    y_mod = (ones(n_mod,1)*y).*exp(i*phase_pro);
    if plotsingle
        figure(2);
        subplot(1,2,1);plot(t,real(y_mod'+ones(length(t),1)*2e-2*(1:n_mod)));
        title('modulated field (real part) in time domain');xlabel('time (ns)');
    end
    y_mod_eff=abs(y_mod);           %this is for the efficiency calculation.
    
    Ys_mod = fftshift(fft(y_mod',n_fourier),1)';
    Ys_mod_eff = fftshift(fft(y_mod_eff',n_fourier),1)';
    Pyy_mod = Ys_mod.* conj(Ys_mod) / n_fourier;
    Pyy_mod_eff = Ys_mod_eff.* conj(Ys_mod_eff) / n_fourier;
    zoom = 10;
    indzoom =round(((1-1/zoom)/2*length(f)):1:((1+1/zoom)/2*length(f)));
    if plotsingle
        subplot(1,2,2);plot(f(indzoom),Pyy_mod(:,indzoom)'+2e-4*ones(length(indzoom),1)*(1:n_mod));
        title('Power spectrum of the modulated pulse');xlabel('frequency (GHz)');
    end
    
    %[param fval error1]= fminsearch(@overlapratio, param,  [],  Pyy_mod,Pyy,f);
    
    %% Calculate the efficiency of detection with a cavity filter.
    err1=sum(Pyy_mod.*(ones(n_mod,1)*Ldistr(f,param)),2)/sum(Pyy);
    eff1=sum(Pyy_mod_eff.*(ones(n_mod,1)*Ldistr(f,param)),2)/sum(Pyy);
    if plotsingle
        figure(3);
        if modtype ==2
            ind = find(num_flip ~= 0);
            plot(num_flip(ind),err1(ind),'bo');
            xlabel('number of flips');
        else if modtype == 3
            ind=2:n_mod;
            plot(ind-1,err1(ind),'bo');
            xlabel('number of shifts');
            end
        end
        title(['PF: eff_{avg}=' num2str(mean(eff1)*100) '%, eff_{min}=' num2str(min(eff1)*100)  '%,   err_{avg}='... 
                num2str(mean(err1(ind))*100) '%, err_{max}=' num2str(max(err1(ind))*100) '%']);
        ylabel('error rate');
    else
        ind = 2:n_mod;
    end
    effavg = [effavg mean(eff1(2:n_mod))];
    effmin = [effmin min(eff1)];
    erravg = [erravg mean(err1(ind))];
    errmax = [errmax max(err1(ind))];
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
    plot(n_mod_range,erravg./effavg,'bo');
    title('BER = Err_{avg}/Eff_{avg}'); xlabel('Encoding Dimension');ylabel('Effective BER');
end

