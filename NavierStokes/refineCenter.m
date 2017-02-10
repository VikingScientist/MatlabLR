function refineCenter(lr, n) 

disp 'Refining geometry';

border = 0.55; % manually chosen, gives resonable mesh around p=2 and 3

for i=1:n
  % collect the bounding box for all basis functions
  xNE = lr.knots(:,lr.p(1)+2);
  yNE = lr.knots(:,end); % north-east
  xNW = lr.knots(:,1);
  yNW = lr.knots(:,end); % nort-west
  xSE = lr.knots(:,lr.p(1)+2);
  ySE = lr.knots(:,lr.p(1)+3); % south-east
  xSW = lr.knots(:,1);
  ySW = lr.knots(:,lr.p(1)+3); % south-west
  g = lr.getGrevillePoint();
  dist = border + sqrt(2);

  % search for all functions around the center
  j = find((xNE-0).^2 + (yNE-0).^2 <= dist^2 & ...
           (xNW-0).^2 + (yNW-0).^2 <= dist^2 & ...
           (xSE-0).^2 + (ySE-0).^2 <= dist^2 & ...
           (xSW-0).^2 + (ySW-0).^2 <= dist^2 );
  j = find(g(:,1).^2 + g(:,2).^2 <= dist^2);
  % search for all functions around the center
  j = j(find((xNE(j)>1 | xNW(j)<-1 | yNE(j)>1 | ySE(j)<-1)));
  
  % manually chosen refinement range decrease
  border = border * 4 / 12;

  % actually do the refinement
  lr.refine(j,'basis');
end
