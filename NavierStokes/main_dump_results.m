
if ~exist(Problem.Title, 'dir')
  mkdir(Problem.Title);
end
filename = sprintf('%s/p%d%d-re%d-T%d', Problem.Title, Problem.Polynomial_Degree, floor(Problem.Reynolds), floor(Problem.Time_Range(2)));
lr.save( [filename, '-lr.lr' ]);
lru.save([filename, '-lru.lr']);
lrv.save([filename, '-lrv.lr']);
lrp.save([filename, '-lrp.lr']);
save(filename, 'Problem', 'uAll', 'time');
