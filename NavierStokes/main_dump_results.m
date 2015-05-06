
if ~exist(Problem.Title, 'dir')
  mkdir(Problem.Title);
end
filename = sprintf('%s/p%d%d-re%d-T%d', Problem.Title, Problem.Polynomial_Degree, floor(Problem.Reynolds), floor(Problem.Time_Range(2)));
save(filename, 'Problem', 'lr', 'lru', 'lrv', 'lrp', 'uAll', 'time');
