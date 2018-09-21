% difference with phasepro: take out the filter part. intended to model the
% mod-demod more precisely. 
function phase_pro = phapro1(dim,t,tau_m,tau_p,method)

if strcmp(method,'pf')
    dim = nextdim(dim);                          %find the next dimision for phase flip.
    phase_pro = zeros(dim, length(t));           %phase profile matrix.
    hada_d=hadamard(dim);
    shuhad = sum(abs(diff(hada_d,1,2)),2)/2;
    [B,Ix]=sort(shuhad);
    %in the order of acsending number of flips.
    phase_mod=hada_d(Ix,:);                            %every possible phase modulation vector.
    %phase_mod=hada_d;

    %Calculate the phase change point for PF. To verify, use int_pul(ind_t).
    int_pul = int_exp(t,tau_p);                  %integrated the lorentzian from -inf to t.
    ind_t = zeros(1,dim);
    for ind_m = 1:dim
        [val,ind_t(ind_m)] = min(abs(int_pul-ind_m/dim));
    end
    phase_pro(:, 1:ind_t(1)) = (phase_mod(:, ones(1,ind_t(1)))-1)/2*pi;
    for ind_m = 1:(dim-1)
        phase_pro(:, ind_t(ind_m):ind_t(ind_m+1)) = (phase_mod(:,(ind_m+1)*ones(1,ind_t(ind_m+1)-ind_t(ind_m)+1))-1)/2*pi;
    end
    
else if strcmp(method,'lp')
        %phase profile matrix.
        phase_pro = pi/tau_m/(dim-1)*(0:(dim-1))'*t;
    else
        error('Unknown modulation type. Use ''pf'' or ''lp''');
    end
end


