function indices = getTensorIndex(lr,el)

i = lr.support{el};
swap_order = [lr.knots(i,lr.p(1)+3:end), lr.knots(i,1:lr.p(1)+2)];
[k, indices] = sortrows(swap_order);
