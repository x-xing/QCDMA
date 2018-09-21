%frequency domain calculation. Used the algebraic results.

% Use qopath to add path.
p1=pwd;
cd([p1(1:strfind(p1, 'Dropbox')-1) 'Dropbox\MatlabCode\utility']);
qopath;
cd(p1);


sup_flag = 1;
if sup_flag 
    a = [0,0.2,0.5,sqrt(0.5)];
else
    a = 0;
end
n_mod = 16;
filt_w = 0.1:0.2:10;

Ixy = zeros(length(a),length(filt_w));
for ind1=1:length(a)
    for ind2 = 1:length(filt_w)
        disp([num2str((ind1-1)*length(filt_w)+ind2) ' out of ' num2str(numel(Ixy))]);
        sp_par=[sup_flag a(ind1)];
        param = [1 n_mod filt_w(ind2)]; %[tau_m, n_mod, filt_w]
        if sp_par(1)
            fid = sp_par(2)^2/(2+sp_par(2)^2); %include the information loss in the symmetric error case.
        else
            fid = 1;
        end
        tic
        % Ixy(ind1,ind2) = cha_cap(param,sp_par)*fid;
        Ixy(ind1,ind2) = cha_cap_asym(param,sp_par); %asymmetric superposition.
        toc
    end
end

figure(sup_flag+1);
pcol = ['bs-';'gx-';'rd-';'k^-'];
hold on;
for ind = 1:4
    plot(filt_w,Ixy(ind,:),pcol(ind,:));
end
hold off;
h=legend('0','0.2','0.5','a');
xlabel('Filter bandwidth (in units of \Delta_{photon})');
ylabel('Mutual information /bits');
title('Mutual information vs. filter');

% Find the first text object (alpha) and change it.
h1 = findobj(get(h,'Children'),'String','a');
set(h1,'String','$\sqrt{2}$','Interpreter','latex')    



