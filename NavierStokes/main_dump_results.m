Problem

if Problem.Static
  if Problem.MatlabPlot
    makePlots;
  end
  
  if Problem.Paraview
    if exist('newElU')==1
      [plotA plotB mesh edges x y] = getSurfacePlotMatrices(lr, lru, lrv, lrp, 4, newElU, newElV, newElP);
    else
      [plotA plotB mesh edges x y] = getSurfacePlotMatrices(lr, lru, lrv, lrp, 4);
    end
    n        = n1+n2;
    u        = uAll;
    vel      = plotA*u(1:n);
    pressure = plotB*u(n+1:end);
    velX     = vel(1:end/2  );
    velY     = vel(  end/2+1:end);
    title = sprintf('%s, %s (%s)', Problem.Title, Problem.Subtitle, Problem.Identifier);
    writeVTK2(strcat(filename, '.vtk'), title,  x,y,mesh,   velX,velY,pressure);
  end
end

if Problem.Save_Results
  filename
  lr.save( [filename, '-lr.lr' ]);
  lru.save([filename, '-lru.lr']);
  lrv.save([filename, '-lrv.lr']);
  lrp.save([filename, '-lrp.lr']);
  save(filename, 'Problem', 'uAll', 'time');
end

