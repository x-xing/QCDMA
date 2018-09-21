function [err indcell]= beropt(Pyx)
dim=size(Pyx,1);

err = zeros(1,dim-1);
indcell = cell(1,dim-1);
for i1 = 2:dim          %encoding dimension
    c1 = nchoosek(1:dim,i1); %all possible sets for this encoding dim.
    er1 = zeros(1,size(c1,1));
    norm = zeros(1,size(c1,1));
    for i2 = 1:size(c1,1) 
        c2 = nchoosek(c1(i2,:),2); %find all possible input-oupt pairs
        c2 = (c2(:,1)-1)*dim+c2(:,2); %translate (row,col) to indices. 
        er1(i2) = mean(Pyx(c2));
        norm(i2) = mean(Pyx((c1(i2,:)-1)*dim+c1(i2,:)));
    end
    er1 = er1./(er1+norm);
    [err(i1-1) ind3] = min(er1);
    indcell{i1-1} = c1(ind3,:);
end

        