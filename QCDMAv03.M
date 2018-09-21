clear all;
close all;
%% Use qopath to add path.
p1=pwd;
cd([p1(1:strfind(p1, 'Dropbox')-1) 'Dropbox\MatlabCode\utility']);
qopath;
cd(p1);

%% Parameter initialization.

t_limit = 200;                      %time axis limit
tau_p = 50;                         %photon lifetime.
tau_m = 10;                         %modulation speed. minimum time to applied Vpi.
%N = 8e3;                            %number of steps 4e4 for PF
N = 20*t_limit/tau_m;               %set the number of steps according to: 10 points per Vpi in LP.
t=linspace(-t_limit, t_limit, N);   %time axis in ns.

param = [1, 1.50, 150];            %[peak transmission, FSR (GHz), finesse] for the filter.
plotsingle = 0;                     %parameter to plot a single dimension, or a range.
modtype = 3;                        %modulation type. 1:UCP, 2: PF, 3: LP
n_mod = 15;                          %dimension of the encoding. valid in single dimension cal.
n_mod_max =64;%round(1/(1-exp(-8*tau_m/tau_p))); %to make sure time bin is larger than eom speed.
n_mod_range = 2:n_mod_max;                 %range of dimension of the encoding. dimension is a variabhhle.
n_fourier = 2^nextpow2(N);                   %dimension of fourier transform.

%typical setting: as a function of dimension: 
%plotsingle = 0; 
%modtype = 3, (linear phase); n_mod_range = 2:15
%modtype = 2, (phase filp); n_mod_range = 2:50.


%% Plot the time and frequency domain of the single photon 
%Lorentzian, exponential in time.
y = exp(-abs(t)/tau_p)/2/tau_p; % field integral is normalized.
%Gaussian
% y = exp(-t.^2/tau_p^2);
figure(1);
subplot(1,2,1);plot(t,y);title('Single photon field in time domain');xlabel('time (ns)')
Ys = fftshift(fft(y,n_fourier));
Pyy = Ys.* conj(Ys) /n_fourier;
f = N/(2*t_limit)*((-n_fourier/2):(n_fourier/2-1))/n_fourier; % see Matlab fft example for the pre-factor.
zoomin = 100;
indzoom =round(((1-1/zoomin)/2*length(f)):((1+1/zoomin)/2*length(f)));
x1 = f(indzoom);
y1 = Pyy(indzoom);
subplot(1,2,2);plot(x1,[y1;fptransm(x1,param)*max(y1)]');  %%%%%%%%%%%%% fptransm power? amplitude????????/
title('Power spectrum of the single photon pulse');xlabel('frequency (GHz)');legend('Photon','Filter');

effavg=[];
effmin=[];
erravg=[];
errmax=[];

if plotsingle
    n_mod_range = n_mod;
end

%% Model the finite EOM speed.
% moving average of the phase profile, with a gaussian kernel.
filt_n = round(tau_m*5*N/t_limit);                     %average of the phase, not the field
filt = 1/(sqrt(2*pi)*filt_n/4)*exp(-((1:filt_n)-round(filt_n/2)).^2/(2*(filt_n/4)^2));
filt=(1-sum(filt))/length(filt)+filt;

%[param fval error1]= fminsearch(@overlapratio, param,  [],
%Pyy_mod,Pyy,f);

%% generate the phase profile
% Phase-flip encoding (PF)
    % Hadamard matrix (n x n) exists only when rem(n,4)=0. 
    % Matlab hadamard(n) function handles only the cases where n, n/12, or n/20 is a power of 2.
if modtype ==2
    d1= min(find(n_mod_max./2.^(1:50)<=1));            %the minimum number to satisfy 2.^d > n_mod.
    d2=100;% min(find(n_mod_max/12./2.^(0:50)<=1));         %the minimum number to satisfy 12*2.^d > n_mod.
    d3= min(find(n_mod_max/20./2.^(0:50)<=1));         %the minimum number to satisfy 20*2.^d > n_mod.
    n_mod_max = min([2^d1; 12*2^(d2-1); 20*2^(d3-1)]); %increase the dimension as required by hadamard.
    if ~plotsingle
        n_mod_range = 2:n_mod_max;
    end
    hada_d=hadamard(n_mod_max);
    phase_mod=hada_d;                                  %every possible phase modulation vector.
    phase_pro = zeros(n_mod_max, length(t));           %phase profile matrix.
    
    %%Calculate the phase change point for PF.
    %To verify, use int_pul(ind_t).
    int_pul = int_exp(t,tau_p);                         %integrated the lorentzian from -inf to t.
    ind_t = zeros(1,n_mod_max);
    for ind_m = 1:n_mod_max
        [val,ind_t(ind_m)] = min(abs(int_pul-ind_m/n_mod_max));
    end
    phase_pro(:, 1:ind_t(1)) = (phase_mod(:, ones(1,ind_t(1)))-1)/2*pi;
    for ind_m = 1:(n_mod_max-1)
        phase_pro(:, ind_t(ind_m):ind_t(ind_m+1)) = (phase_mod(:,(ind_m+1)*ones(1,ind_t(ind_m+1)-ind_t(ind_m)+1))-1)/2*pi;
    end
    
    phase_pro = [phase_pro(:,1)*ones(1,filt_n) phase_pro];  %padding the first filt_n datapoints of phase_pro.
    phase_pro=filter(filt,1,phase_pro,[],2);                %moving average filtering.
    phase_pro(:,1:filt_n)=[];                               %remove the padding.

    y_mod = (ones(n_mod_max,1)*y).*exp(i*phase_pro);
    y_mod_eff=abs(y_mod);           %this is for the efficiency calculation.
    
    Ys_mod = fftshift(fft(y_mod',n_fourier),1)';
    Ys_mod_eff = fftshift(fft(y_mod_eff',n_fourier),1)';
    Pyy_mod = Ys_mod.* conj(Ys_mod) / n_fourier;
    Pyy_mod_eff = Ys_mod_eff.* conj(Ys_mod_eff) / n_fourier;
    
    err1=sum(Pyy_mod.*(ones(n_mod_max,1)*fptransm(f,param)),2)/sum(Pyy);
    eff1=sum(Pyy_mod_eff.*(ones(n_mod_max,1)*fptransm(f,param)),2)/sum(Pyy);
    
    %pick the minimum error rate encoding.
    [err1_st, ind_st]=sort(err1,'ascend');
%     %%%in order of descending number of flips
%      p1 = diff(phase_mod,1,2); %difference of 1 along dimension 2 (coloumn).
%      [nf_st, ind_st]=sort(sum(abs(p1)/2,2),'descend'); %sort according to number of flips.
%    ind_st(2) = 2;

end

for n_mod =n_mod_range
    sprintf('%d of %d',n_mod,max(n_mod_range))
    %1. Unit-circle phase encoding (UCP)
    %to-do
    if modtype ==1                 
        phase_mod = exp(i*2*pi*(1:n_mod)/n_mod);           %phase changes in steps.
        error('Unit circle phase encoding is not implemented yet!');
    end
        
    %3. Linear phase encoding (LP)
    if modtype ==3
        phase_pro = zeros(n_mod, length(t));         %phase profile matrix.
        phase_pro = pi/tau_m/(n_mod-1)*(0:(n_mod-1))'*t;
        phase_pro = [phase_pro(:,1)*ones(1,filt_n) phase_pro];        %padding the first filt_n datapoints of phase_pro.
        phase_pro = filter(filt,1,phase_pro,[],2);                      %moving average filtering.
        phase_pro(:,1:filt_n) = [];                                     %remove the padding.
        y_mod = (ones(n_mod,1)*y).*exp(i*phase_pro);
        y_mod_eff = abs(y_mod);           %this is for the efficiency calculation.

        Ys_mod = fftshift(fft(y_mod',n_fourier),1)';
        Ys_mod_eff = fftshift(fft(y_mod_eff',n_fourier),1)';
        Pyy_mod = Ys_mod.* conj(Ys_mod) / n_fourier;
        Pyy_mod_eff = Ys_mod_eff.* conj(Ys_mod_eff) / n_fourier;

        %Calculate the efficiency of detection with a cavity filter.
        err1 = sum(Pyy_mod.*(ones(n_mod,1)*fptransm(f,param)),2)/sum(Pyy);
        eff1 = sum(Pyy_mod_eff.*(ones(n_mod,1)*fptransm(f,param)),2)/sum(Pyy);

        effavg = [effavg mean(eff1)];
        effmin = [effmin min(eff1)];
        erravg = [erravg mean(err1(2:n_mod))];
        errmax = [errmax max(err1(2:n_mod))];        
    end
    
    if plotsingle
        figure(2);
        if modtype ==2
%             y_mod = y_mod(ind_st(1:4:n_mod),:); for the QCDMA paper
%             Pyy_mod = Pyy_mod(ind_st(1:4:n_mod),:);
             y_mod = y_mod(ind_st(1:(n_mod-1)),:);
             Pyy_mod = Pyy_mod(ind_st(1:(n_mod-1)),:);
            subplot(1,2,1);plot(t,real(y_mod'+ones(length(t),1)*4e-2*(1:length(y_mod(:,1)))),t,real(y));
            title('modulated field (real part) in time domain');xlabel('time (ns)');
            zoom = 100;
            indzoom =round(((1-1/zoom)/2*length(f)):1:((1+1/zoom)/2*length(f)));
            subplot(1,2,2);plot(f(indzoom),Pyy_mod(:,indzoom)'+2.5e-2*ones(length(indzoom),1)*(1:length(y_mod(:,1))),f(indzoom),Pyy(:,indzoom)'/6);
            hold on;
            plot(f(indzoom),max(Pyy_mod(1,indzoom))*fptransm(f(indzoom),param),'k-.');hold off;
            title('Power spectrum of the modulated pulse');xlabel('frequency (GHz)');
        else if modtype == 3
                subplot(1,2,1);plot(t,real(y_mod'+ones(length(t),1)*2e-2*(1:length(y_mod(:,1)))));
                title('modulated field (real part) in time domain');xlabel('time (ns)');
                zoom = 100;
                indzoom =round(((1-1/zoom)/2*length(f)):1:((1+1/zoom)/2*length(f)));
                subplot(1,2,2);plot(f(indzoom),Pyy_mod(:,indzoom)'+2.5e-2*ones(length(indzoom),1)*(0:(length(y_mod(:,1))-1)));
                hold on;
                plot(f(indzoom),max(Pyy_mod(1,indzoom))*fptransm(f(indzoom),param),'k-.');hold off;
                title('Power spectrum of the modulated pulse');xlabel('frequency (GHz)');
            end
        end
        

        figure(3);
        if modtype ==2
            if n_mod == n_mod_max          % not to consider 0 flips in error cal.
                ind = ind_st(1:(n_mod-1));
            else
                ind = ind_st(1:n_mod);
            end
            plot(err1(ind),'bo');
            xlabel('index');
            title(['PF: eff_{avg}=' num2str(mean(eff1)*100) '%, eff_{min}=' num2str(min(eff1)*100)  '%,   err_{avg}='... 
                num2str(mean(err1(ind))*100) '%, err_{max}=' num2str(max(err1(ind))*100) '%']);
        else if modtype == 3
            ind=2:n_mod;
            plot(ind-1,err1(ind),'bo');
            xlabel('number of shifts');
            title(['PF: eff_{avg}=' num2str(mean(eff1)*100) '%, eff_{min}=' num2str(min(eff1)*100)  '%,   err_{avg}='... 
                num2str(mean(err1(ind))*100) '%, err_{max}=' num2str(max(err1(ind))*100) '%']);
            end
        end
        ylabel('error rate');
    end
    
    if modtype == 2
        if n_mod == n_mod_max          % not to consider 0 flips in error cal.
            ind = ind_st(1:(n_mod-1));
        else
            ind = ind_st(1:n_mod);
        end
        erravg = [erravg mean(err1(ind))];
        errmax = [errmax max(err1(ind))];
        effavg = [effavg mean(eff1(ind))];
        effmin = [effmin min(eff1(ind))];
    end
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
