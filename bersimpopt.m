function err = bersimpopt(Pyx,ratio_tm_tp)
dim= size(Pyx,1);
err = zeros(1,dim-1);
for d = 2:dim
    %v=round(linspace(1,dim,d));
    %maximize the XOR value. 
    v=vtoind(vopt(dim,d),ratio_tm_tp);
    %check xor value of v
    %sumxor(v)
    %Was the following: incorrect.
%     c = nchoosek(v,2); %find every possible combination of input-output.
%     c = (c(:,1)-1)*dim+c(:,2); %translate (row,col) to absolute indices
%     er1 = mean(Pyx(c));%average error probability
%     norm = mean(Pyx((v-1)*dim+v)); %normalized to diagonal probabilities
%     err(d-1)= er1./(er1+norm);
    for i1 = 1:d
        pii = Pyx(v(i1),v(i1));
        sumpij = sum(Pyx(v(i1),v));
        err(d-1)=err(d-1)+1/d*(1-pii/sumpij);
    end
end