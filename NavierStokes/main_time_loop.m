
nSteps = length(time);
% time = linspace(0,30,nSteps);
% k = time(2)-time(1);
k = delta_time;
N = n1+n2+n3;
u        = zeros(N,1);
u(velEdges)    = velVal;
u(presEdges+n) = presVal;
% u(lru.getEdge(4)) = 1;
uAll = zeros(N,nSteps);
uAll(:,1) = u;


% lhs = [M+k*A, k*D; D', zeros(n3,n3)];
% lhs = [M+k*A, k*D; Dt, [avg_p'; zeros(n3-1,n3)]];
% lhs = [M, zeros(n1+n2,n3); zeros(n3,N)] + k*dF(1);
% b   = b(1:(n1+n2));
% lhs(edges,:) = [];
% lhs(:,edges) = [];

% [plotAu meshu eu xu yu] = lru.getSurfMatrix('diffX', 'parametric', 'nviz', 5, 'diffX');
% [plotAv meshu ev xv yv] = lrv.getSurfMatrix('diffY', 'parametric', 'nviz', 5, 'diffY');
if exist('newElU')==1
  [plotA plotB mesh edges x y] = getSurfacePlotMatrices(lr, lru, lrv, lrp, 4, newElU, newElV, newElP);
else
  [plotA plotB mesh edges x y] = getSurfacePlotMatrices(lr, lru, lrv, lrp, 4);
end
title    = sprintf('%s, %s (%s)', Problem.Title, Problem.Subtitle, Problem.Identifier);
writeVTK2(sprintf('%s-%d.vtk', filename, 0), title, x,y,mesh, zeros(numel(x),1),zeros(numel(x),1), zeros(numel(x),1));

topRightCorner = intersect(lrp.getEdge(2), lrp.getEdge(4))+n1+n2;

% b = -A(:,velEdges)*velVal - D(:,presEdges)*presVal;
timer = cputime; tic;
for i=2:nSteps
  lastTime = toc;
  fprintf('Time: %g (step %d/%d):\n', time(i), i, nSteps);

  % initial guess for newton stepping = previous time step
  v   = u; %zeros(n1+n2+n3,1);
  n  = n1+n2;
  N  = n1+n2+n3;
  % rhs = M*(v(1:n)-u(1:n)) + k*(A*v(1:n)+D*v(n+1:end));
  % rhs = M*u(1:n) - k*b - M(:,velEdges)*velVal;
  % lhs = [M + k*A, k*D];
  % lhs(velEdges,:) = 0;
  % lhs(:,velEdges) = 0;
  % lhs(velEdges,velEdges) = speye(numel(velEdges));

  % rhs(velEdges) = velVal;
  lhs_lower = [D', zeros(n3)];
  lhs_lower([presEdges;avgP_ind],:) = 0;
  lhs_lower(avgP_ind,n+1:end)      = avg_p';
  lhs_lower(presEdges,n+presEdges) = speye(numel(presEdges));;
  % rhs_lower = bLow;
  rhs_lower = zeros(n3,1);
  rhs_lower(presEdges) = presVal;
  rhs_lower(avgP_ind)  = 0;
  % rhs = [rhs; -D(velEdges,:)'*velVal];
  % rhs(presEdges+n) = presVal;
  % rhs(avgP_ind +n) = 0;
  lhs = [[M,zeros(n,n3)] + k*dF(zeros(size(u))); lhs_lower];
  rhs = [M*u(1:n)        + k* F(zeros(size(u))); rhs_lower];
  v = lhs \ rhs;
  for newtIt=1:nwtn_max_it
    %%% QUACK!!!
    rhs = M*(v(1:n)-u(1:n)) + k*(A*v(1:n)+D*v(n+1:end));
    % rhs = rhs - A(:,velEdges)*velVal - D(:,presEdges)*presVal;
    lhs = [M + k*A, k*D];
    lhs(velEdges,:) = 0;
    lhs(:,velEdges) = 0;
    % for j=1:numel(velEdges)
      % k = velEdges(j);
      % lhs(k,k) = 1;
    % end
    lhs(velEdges,velEdges) = speye(numel(velEdges));

    rhs(velEdges) = 0;
    lhs_lower = [D', zeros(n3)];
    lhs_lower([presEdges;avgP_ind],:) = 0;
    lhs_lower(avgP_ind,n+1:end)      = avg_p';
    lhs_lower(presEdges,n+presEdges) = speye(numel(presEdges));;
    rhs = [rhs; -D(velEdges,:)'*velVal];
    rhs(presEdges+n) = -presVal;
    rhs(avgP_ind +n) = 0;
    lhs = [lhs; lhs_lower];
    %%% backward euler stepping
    % lhs = [M, zeros(n,n3);    zeros(n3,N)] + k*dF(v);
    % rhs = [M*(v(1:n)-u(1:n)); zeros(n3,1)] + k* F(v);
    %%% crank-nicolson rule stepping
    % lhs = [M, zeros(n,n3);    zeros(n3,N)] + k/2*(dF(v)       );
    % rhs = [M*(v(1:n)-u(1:n)); zeros(n3,1)] + k/2*( F(v)+ F(u) );
    % max(max(abs(rightLHS-lhs)))
    % max(max(abs(rightRHS-rhs)))
    % rhs(edges) = 0;
    dv = lhs \ -rhs;
    v(1:n) = v(1:n) + dv(1:n);
    v(n+1:end) = v(n+1:end);
    if(norm(dv(1:n))<nwtn_res_tol)
      break;
    end
  end
  fprintf('  Newton iteration converged after %d iterations at residual %g\n', newtIt, norm(dv));
  fprintf('  Pressure at top right corner : %g\n', v(topRightCorner));
  fprintf('  Max u-velocity controlpoint  : %g\n', max(v(1:n1)));
  fprintf('  Max v-velocity controlpoint  : %g\n', max(v(n1+1:n1+n2)));
  fprintf('  Max pressure controlpoint    : %g\n', max(v(n+1:end)));
  u = v;
  vel = plotA*u(1:n);
  pressure = plotB*u(n+1:end);
  velX     = vel(1:end/2  );
  velY     = vel(  end/2+1:end);
  writeVTK2(sprintf('%s-%d.vtk', filename, i-1), title, x,y,mesh, velX,velY,pressure);


  % u = [M+k/2*A, k/2*D; Dt [avg_p'; zeros(n3-1,n3)]] \ [(M-k/2*A)*u(1:n)-k/2*D*u(n+1:end)+k*b; zeros(n3,1)];
  % u = [M+k*A, k*D; Dt [avg_p'; zeros(n3-1,n3)]] \ [M*u(1:n)+k*b; zeros(n3,1)];
  % u = rightLHS \ rightRHS;
  % dudx = plotAu*u(1:n1);
  % dvdy = plotAv*u(n1+1:n1+n2);
  % divPt = dudx + dvdy;
  % fprintf('  max(div(u)) = %g\n', max(max(divPt)));
  timeUsed = toc - lastTime;
  timeLeft = (nSteps-i)*timeUsed;
  c = clock;
  c(6) = c(6) + floor(timeLeft);
  fprintf('  Estimated time left: ');
  if timeLeft>60*60*24
    fprintf('%d days, ', floor(timeLeft/60/60/24));
    timeLeft = mod(timeLeft,60*60*24);
  end
  if timeLeft>60*60
    fprintf('%d hours, ', floor(timeLeft/60/60));
    timeLeft = mod(timeLeft,60*60);
  end
  if timeLeft>60
    fprintf('%d minutes and ', floor(timeLeft/60));
    timeLeft = mod(timeLeft,60);
  end
  fprintf('%d seconds\n', floor(timeLeft));
  fprintf('  Estimated completion: %s\n', datestr(datetime(c)));

  uAll(:,i) = u;
end
fprintf('\n');
time_timeStepping = cputime - timer; time_timeStepping_wall = toc;

