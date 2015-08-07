Problem

if Problem.MatlabPlot
  makePlots;
end

if Problem.Save_Results
  filename
  lr.save( [filename, '-lr.lr' ]);
  lru.save([filename, '-lru.lr']);
  lrv.save([filename, '-lrv.lr']);
  lrp.save([filename, '-lrp.lr']);
  save(filename, 'Problem', 'uAll', 'time');
end

