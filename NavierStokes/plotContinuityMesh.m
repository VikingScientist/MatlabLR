function H = plotContinuityMesh(lr)
% function H = plotContinuityMesh(lr)

H = gcf;

n = size(lr.lines,1);
p = lr.p;
TOL = 1e-10;
i = find(abs(lr.lines(:,1) - lr.lines(:,3))<TOL); % const u lines (vertical)
j = find(abs(lr.lines(:,2) - lr.lines(:,4))<TOL); % const v lines (horizontal)

% continuity_colors = [0,1,1; 1/8,7/8,1; 2/8,6/8,1; .6,5/8,.8; 6/8,6/8,1; 1,0,7/8];
continuity_colors = [0,1,1; 2/8,6/8,1; .6,5/8,.8; 1/8,7/8,1; 6/8,6/8,1; 1,0,7/8];

holdstate = ishold;
hold on;
tags = cell(max(p)+1,1);
for c=-1:max(p)-1, % loop over all possible continuities
  if c==-1,
    plot(0,0, 'k-');
  else
    plot(0,0, 'Color', hsv2rgb(continuity_colors(c+1,:)));
  end
  tags{c+2} = sprintf('C^{%d}-lines', c);
end
legend(tags);


for c=p(1)-1:-1:-1, % loop over all vertical continuities
  k = find(c == p(1) - lr.lines(i,5));
  x = lr.lines(i(k), [1,3]);
  y = lr.lines(i(k), [2,4]);
  if c==-1,
    plot(x', y', 'k-');
  else
    plot(x', y', 'Color', hsv2rgb(continuity_colors(c+1,:)));
  end
end
for c=p(2)-1:-1:-1, % loop over all horizontal continuities
  k = find(c == p(2) - lr.lines(j,5));
  x = lr.lines(j(k), [1,3]);
  y = lr.lines(j(k), [2,4]);
  if c==-1,
    plot(x', y', 'k-');
  else
    plot(x', y', 'Color', hsv2rgb(continuity_colors(c+1,:)));
  end
end
elMin = min(lr.elements);
elMax = max(lr.elements);
du    = elMax(3:4) - elMin(1:2);
% axis([elMin(1)-du(1)*.1, elMax(3)+du(1)*.1, elMin(2)-du(2)*.1, elMax(4)+du(2)*.1]);
% axis equal;

if ~holdstate
  hold off;
end


