if iteration_h > 1
  lr.refine();
  h_val_result = h_val_result / 2;
  disp 'Deriving spaces'
  [lrv lru] = lr.getDerivative( 'no cp');
  [lrp   ~] = lru.getDerivative('no cp');
  if doCrop
    [oldBasis  oldEl ] = lr.clipArea(crop);
    [oldBasisU oldElU] = lru.clipArea(crop);
    [oldBasisV oldElV] = lrv.clipArea(crop);
    [oldBasisP oldElP] = lrp.clipArea(crop);
    newElU = zeros(size(lru.elements,1),1);
    newElV = zeros(size(lrv.elements,1),1);
    newElP = zeros(size(lrp.elements,1),1);
    newElU(oldElU) = 1:numel(oldElU);
    newElV(oldElV) = 1:numel(oldElV);
    newElP(oldElP) = 1:numel(oldElP);
  end
end

