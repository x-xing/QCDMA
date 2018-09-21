%generate a movie to display the choice of basis

check = exist('Pyx','var') && exist('indcell','var');
if check
    dim = size(Pyx,1);
    maxP = max(reshape(Pyx,numel(Pyx),[]));
    for i1 = 1:dim-1
        Pyx1=Pyx;
        pcolor(Pyx);colormap('hot');colorbar;
        mov1(2*i1-1)=getframe;
        c=nchoosek(indcell{i1},2);
        c = (c(:,1)-1)*dim+c(:,2);
        Pyx1(c)=maxP*ones(length(c),1);
        pcolor(Pyx1);colormap('hot');colorbar;
        mov1(2*i1)=getframe;
    end
end
figure(4);
movie(mov1,3,1);