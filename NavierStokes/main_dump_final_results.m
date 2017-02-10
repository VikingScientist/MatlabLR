Problem
u = uAll;

t = cputime; tic;
if exist('Convergence_rates')

  yStart = result_uh_H1(end, :);
  yStart = yStart * 1.3;
  xStart = result_h(end,1);
  xRange = 2^(Convergence_rates.iterations-1);
  expected_rates = Convergence_rates.p_values;
  helpLinesY = [yStart ; yStart .*xRange.^expected_rates];
  helpLinesX = [xStart ; xStart  *xRange];
  legendEntries = cell(0);
  for i=1:numel(expected_rates)
    legendEntries{i} = sprintf('p=%d', Convergence_rates.p_values(i));
  end
  for i=1:numel(expected_rates)
    legendEntries{i+numel(expected_rates)} = sprintf('$$\\mathcal{O}(h^{-%d})$$', expected_rates(i));
  end
  figure;
    loglog(result_h, result_uh_H1, 'o-');
    c = get(gca, 'ColorOrder');
    set(    gca, 'ColorOrder', c(1:numel(expected_rates),:));
    hold on;
    loglog(helpLinesX, helpLinesY, '--');
    title('Velocity error');
    legend(legendEntries, 'interpreter', 'latex', 'Location', 'NorthWest');
    set(gca, 'FontSize', 16);
    xlabel('Mesh size $$h$$', 'interpreter', 'latex')
    ylabel('$$\|u-u_h\|_{H^1} / \|u\|_{H^1}$$', 'interpreter', 'Latex');
    % saveas(gcf, 'conv-umr-u-cavity-re200000-cheat.pdf', 'pdf');

  yStart = result_ph_L2(end, :);
  yStart = yStart * 1.5;
  xStart = result_h(end,1);
  xRange = 2^(Convergence_rates.iterations-1);
  expected_rates = Convergence_rates.p_values + 1;
  helpLinesY = [yStart ; yStart .*xRange.^expected_rates];
  helpLinesX = [xStart ; xStart  *xRange];
  for i=1:numel(expected_rates)
    legendEntries{i} = sprintf('p=%d', Convergence_rates.p_values(i));
  end
  for i=1:numel(expected_rates)
    legendEntries{i+numel(expected_rates)} = sprintf('$$\\mathcal{O}(h^{-%d})$$', expected_rates(i));
  end
  figure;
    loglog(result_h, result_ph_L2, 'o-');
    c = get(gca, 'ColorOrder');
    set(    gca, 'ColorOrder', c(1:numel(expected_rates),:));
    hold on;
    loglog(helpLinesX, helpLinesY, '--');
    title('Pressure error');
    legend(legendEntries, 'interpreter', 'latex', 'Location', 'NorthWest');
    set(gca, 'FontSize', 16);
    xlabel('Mesh size $$h$$', 'interpreter', 'latex')
    ylabel('$$\|p-p_h\|_{L^2} / \|p\|_{L^2}$$', 'interpreter', 'Latex');
    % saveas(gcf, 'conv-umr-p-cavity-re200000-cheat.pdf', 'pdf');
end

if Problem.Static
  if Problem.MatlabPlot
    % makePlots;

    n        = n1+n2;
    u        = uAll;
    if exist('newElU')==1
      [plotA plotB mesh edges x y] = getSurfacePlotMatrices(lr, lru, lrv, lrp, 4, newElU, newElV, newElP);
    else
      [plotA plotB mesh edges x y] = getSurfacePlotMatrices(lr, lru, lrv, lrp, 4);
    end
    vel      = plotA*u(1:n);
    pressure = plotB*u(n+1:end);
    velX     = vel(1:end/2  );
    velY     = vel(  end/2+1:end);
    xrange   = [min(x), max(x)];
    yrange   = [min(y), max(y)];

    z = velX;
    figure; hold on;
      patch('Faces', mesh, 'Vertices', [x,y,z], 'CData', z, 'FaceColor', 'interp', 'EdgeColor', 'none');
      plot3(x(edges), y(edges), z(edges), 'k-');
      set(gca, 'FontSize', 24);
      title('U solution');
      xlim(xrange);
      ylim(yrange);
      colorbar;
      axis equal;
    z = velY;
    figure; hold on;
      patch('Faces', mesh, 'Vertices', [x,y,z], 'CData', z, 'FaceColor', 'interp', 'EdgeColor', 'none');
      plot3(x(edges), y(edges), z(edges), 'k-');
      set(gca, 'FontSize', 24);
      title('V solution');
      xlim(xrange);
      ylim(yrange);
      colorbar;
      axis equal;
    z = pressure;
    figure; hold on;
      patch('Faces', mesh, 'Vertices', [x,y,z], 'CData', z, 'FaceColor', 'interp', 'EdgeColor', 'none');
      plot3(x(edges), y(edges), z(edges), 'k-');
      set(gca, 'FontSize', 24);
      title('P solution');
      xlim(xrange);
      ylim(yrange);
      colorbar;
      axis equal;
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
time.plot     = cputime - t;
walltime.plot = toc;

t = cputime; tic;
if Problem.Save_Results
  filename
  lr.save( [filename, '-lr.lr' ]);
  lru.save([filename, '-lru.lr']);
  lrv.save([filename, '-lrv.lr']);
  lrp.save([filename, '-lrp.lr']);
  if Problem.Static
    save(filename, 'Problem', 'uAll');
  else
    save(filename, 'Problem', 'uAll', 'time');
  end
end
time.savetofile     = cputime - t;
walltime.savetofile = toc;

