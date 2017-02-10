function refineCorners(lr, n) 

disp 'refining LR spline';
dist = 5/16;
hmin = 1/16;
for k=1:n,
  dist = hmin*6.5;
  x = lr.knots(:,lr.p(1)+2);
  y = lr.knots(:,end);
  i = find( (x+0).^2 + (y+0).^2 <  dist^2);
  x = lr.knots(:,1);
  y = lr.knots(:,end);
  j = find( (x-1).^2 + (y+0).^2 <  dist^2);
  x = lr.knots(:,lr.p(1)+2);
  y = lr.knots(:,lr.p(1)+3);
  ii = find( (x+0).^2 + (y-1).^2 < dist^2);
  x = lr.knots(:,1);
  y = lr.knots(:,lr.p(1)+3);
  jj = find( (x-1).^2 + (y-1).^2 < dist^2);
  lr.refine([i;j;ii;jj], 'basis');
  dist = 2*dist / 3;
  hmin = hmin/2;
  % figure; lr.plot('parametric'); axis equal;
  % saveas(gcf, sprintf('Results/mesh-%d.png', k), 'png');
  % saveas(gcf, sprintf('Results/mesh-%d.eps', k), 'psc2');
end

