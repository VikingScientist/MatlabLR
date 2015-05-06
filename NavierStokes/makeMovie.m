function makeMovie(filename)

lr = LRSplineSurface
lru = LRSplineSurface
lrv = LRSplineSurface
lrp = LRSplineSurface

load(filename);
lr.load( [filename, '-lr.lr' ]);
lru.load([filename, '-lru.lr']);
lrv.load([filename, '-lrv.lr']);
lrp.load([filename, '-lrp.lr']);

nSteps = length(time);

% pre-compute all matrices needed to plot the solutions
[plotA mesh edges x y]      = getSurfacePlotMatrices(lr, lru, lrv, lrp, 5);
[Aquiv xquiv, yquiv, Jquiv] = getQuiverPlotMatrices(lru, lrv, 30,30, lr);

% initialize movie capture
myMovie(nSteps) = struct('cdata', [], 'colormap', []);
figure;
set(gcf, 'Position', [0,0, 1280, 800]);

for i=1:nSteps
  u = uAll(1:n1,i);
  v = uAll((n1+1):(n1+n2),i);
  vel      = Aquiv*[u;v];
  quivVelX = vel(1:end/2);
  quivVelY = vel(end/2+1:end);
  vel  = plotA*[u;v];
  velX = vel(1:end/2);
  velY = vel(end/2+1:end);
  z    = sqrt(velX.^2 + velY.^2);
  clf; hold on;
    patch('Faces', mesh, 'Vertices', [x,y,zeros(size(x))], 'CData', z, 'FaceColor', 'interp', 'EdgeColor', 'none');
    quiver(xquiv, yquiv, quivVelX, quivVelY, 4.0, 'LineWidth', 2, 'Color', 'Black');
    plot3(x(edges), y(edges), zeros(size(edges)), 'k-');
    colorbar;
    set(gca, 'CLim', [-.3, 1.0]);
    xlim([0 1]);
    ylim([0 1]);
    title(sprintf('Time t=%.3f', time(i)), 'FontSize', 24);
  myMovie(i) = getframe;
end

movie(myMovie);
movie2avi(myMovie, [filename, '.avi']);

