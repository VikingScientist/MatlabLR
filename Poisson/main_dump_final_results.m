Problem
u = uAll;

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
    legendEntries{i+numel(expected_rates)} = sprintf('$$\\mathcal{O}(h^{%d})$$', expected_rates(i));
  end
  figure;
    loglog(result_h, result_uh_H1, 'o-');
    c = get(gca, 'ColorOrder');
    set(    gca, 'ColorOrder', c(1:numel(expected_rates),:));
    hold on;
    loglog(helpLinesX, helpLinesY, '--');
    title('Energy error');
    legend(legendEntries, 'interpreter', 'latex', 'Location', 'NorthWest');
    set(gca, 'FontSize', 16);
    xlabel('Mesh size $$h$$', 'interpreter', 'latex')
    ylabel('$$|u-u_h|_{H^1} / |u|_{H^1}$$', 'interpreter', 'Latex');
    % saveas(gcf, 'buffa-conv-u-zero-p-in-corner.pdf', 'pdf');

  yStart = result_uh_L2(end, :);
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
    legendEntries{i+numel(expected_rates)} = sprintf('$$\\mathcal{O}(h^{%d})$$', expected_rates(i));
  end
  figure;
    loglog(result_h, result_uh_L2, 'o-');
    c = get(gca, 'ColorOrder');
    set(    gca, 'ColorOrder', c(1:numel(expected_rates),:));
    hold on;
    loglog(helpLinesX, helpLinesY, '--');
    title('L2 error');
    legend(legendEntries, 'interpreter', 'latex', 'Location', 'NorthWest');
    set(gca, 'FontSize', 16);
    xlabel('Mesh size $$h$$', 'interpreter', 'latex')
    ylabel('$$\|u-u_h\|_{L^2} / \|u\|_{L^2}$$', 'interpreter', 'Latex');
    % saveas(gcf, 'buffa-conv-p-zero-p-in-corner.pdf', 'pdf');
end

if Problem.Static
  if Problem.MatlabPlot
    figure; lr.surf(u);         title('u');;
    % figure; lr.surf(u, 'diffX');title('du/dx');
    % figure; lr.surf(u, 'diffY');title('du/dy');
    figure; lr.surf(u, 'secondary', @(x,u,du) norm(du-Exact_solution.grad_u(x(1),x(2))), 'nviz', 10); title('Enrgy err');
    view(2); colorbar;
    figure; lr.surf(u, 'secondary', @(x,u,du) norm(u-Exact_solution.u(x(1),x(2))), 'nviz', 10); title('L2 error');
    view(2); colorbar;
    figure; lr.surf(lr.p(:,1)); title('polynomial degree'); colorbar;

  end
end

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

