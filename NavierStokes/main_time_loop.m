
nSteps = length(time);
k = delta_time;
N = n1+n2+n3;
n = n1+n2;
u        = zeros(N,1);
edg = [velEdges; presEdges+n];
val = [velVal  ; presVal];
u(edg) = val*timeScale(0);
nonEdge = 1:N;
nonEdge(edg) = [];

uAll = zeros(N,nSteps);
uAll(:,1) = u;

ndof_vel = n1+n2-numel(velEdges);


if Problem.Paraview
  % [plotAu meshu eu xu yu] = lru.getSurfMatrix('diffX', 'parametric', 'nviz', 5, 'diffX');
  % [plotAv meshu ev xv yv] = lrv.getSurfMatrix('diffY', 'parametric', 'nviz', 5, 'diffY');
  if exist('newElU')==1
    [plotA plotB mesh edges x y] = getSurfacePlotMatrices(lr, lru, lrv, lrp, 4, newElU, newElV, newElP);
  else
    [plotA plotB mesh edges x y] = getSurfacePlotMatrices(lr, lru, lrv, lrp, 4);
  end
  title    = sprintf('%s, %s (%s)', Problem.Title, Problem.Subtitle, Problem.Identifier);
  writeVTK2(sprintf('%s-%d.vtk', filename, 0), title, x,y,mesh, zeros(numel(x),1),zeros(numel(x),1), zeros(numel(x),1));
end

topRightCorner = intersect(lrp.getEdge(2), lrp.getEdge(4))+n1+n2;

timer = cputime; tic;
for i=2:nSteps
  lastTime = toc;
  fprintf('Time: %g (step %d/%d):\n', time(i), i, nSteps);
  integratorUsed = '';
  tnp1 = time(i);   % t_{n+1}
  tn   = time(i-1); % t_{n}
  dt   = tnp1-tn;

  % initial guess for newton stepping = previous time step
  v   = u;

  for newtIt=1:Problem.Newton_Max_It
    un = u(nonEdge);  % u_{n}
    vn = v(nonEdge);  % u_{n+1} (newton-iteration approximation of it)

    if i<0 || strcmp(time_integrator, 'Backward Euler')
      integratorUsed = 'Backward Euler';
      lhs = dt*dF(vn,tnp1);
      rhs = dt* F(vn,tnp1);
    elseif strcmp(time_integrator, 'Crank-Nicolsen')
      integratorUsed = 'Crank-Nicolsen';
      lhs = dt/2*(dF(vn,tnp1)            );
      rhs = dt/2*( F(vn,tnp1) + F(un,tn) );
    end

    lhs(1:ndof_vel,1:ndof_vel) = lhs(1:ndof_vel,1:ndof_vel) + M;
    rhs(1:ndof_vel)            = rhs(1:ndof_vel)            + M*(vn(1:ndof_vel)-un(1:ndof_vel));
    dv = lhs \ -rhs;

    v(nonEdge) = v(nonEdge) + dv;
    if(norm(dv)<Problem.Newton_TOL)
      break;
    end
  end
  fprintf('  Newton iteration converged after %d iterations at residual %g\n', newtIt, norm(dv));
  fprintf('  Time integrator              : %s\n', integratorUsed);
  fprintf('  Pressure at top right corner : %g\n', v(topRightCorner));
  fprintf('  Max u-velocity controlpoint  : %g\n', max(v(1:n1)));
  fprintf('  Max v-velocity controlpoint  : %g\n', max(v(n1+1:n1+n2)));
  fprintf('  Max pressure controlpoint    : %g\n', max(v(n+1:end)));
  u = v;
  u(edg) = val*timeScale(tnp1);
  if Problem.Paraview
    vel = plotA*u(1:n);
    pressure = plotB*u(n+1:end);
    velX     = vel(1:end/2  );
    velY     = vel(  end/2+1:end);
    writeVTK2(sprintf('%s-%d.vtk', filename, i-1), title, x,y,mesh, velX,velY,pressure);
  end

  % compute and print timing estimates
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
  % fprintf('  Estimated completion: %s\n', datestr(datetime(c)));

  uAll(:,i) = u;
end
fprintf('\n');
time_timeStepping = cputime - timer; time_timeStepping_wall = toc;

