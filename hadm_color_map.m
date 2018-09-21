ha16=hadamard(16);
[hs ind]=sort(sum(abs(diff(ha16,1,2)),2),'ascend');

pcolor(ha16(ind,:));
colormap(gray(2));
xlabel('Column');
ylabel('Row');
title('Colour map of the Walsh matrix');