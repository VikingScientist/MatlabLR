

n1 = size(lru.knots,1);
n2 = size(lrv.knots,1);
n3 = size(lrp.knots,1);

%%% preassure field boundary conditions
e1 = lrp.getEdge(1);
e2 = lrp.getEdge(2);
e3 = lrp.getEdge(3);
e4 = lrp.getEdge(4);
c1 = intersect(e1,e3);
c2 = intersect(e1,e4);
c3 = intersect(e2,e3);
c4 = intersect(e2,e4);

% c1 = find(lrp.knots(:,2)==-1);
% c1 = c1(find(lrp.knots(c1,6)==-3));
% c1 = c1(find(lrp.knots(c1,3)~=-1));
% c1 = c1(find(lrp.knots(c1,7)~=-3));
% 
% c2 = find(lrp.knots(:,3)== 1);
% c2 = c2(find(lrp.knots(c2,6)==-3));
% c2 = c2(find(lrp.knots(c2,2)~= 1));
% c2 = c2(find(lrp.knots(c2,7)~=-3));
% 
% c3 = find(lrp.knots(:,3)== 1);
% c3 = c3(find(lrp.knots(c3,7)== 3));
% c3 = c3(find(lrp.knots(c3,2)~= 1));
% c3 = c3(find(lrp.knots(c3,6)~= 3));
% 
% c4 = find(lrp.knots(:,2)==-1);
% c4 = c4(find(lrp.knots(c4,7)== 3));
% c4 = c4(find(lrp.knots(c4,3)~=-1));
% c4 = c4(find(lrp.knots(c4,6)~= 3));

corners = [c1,c2,c3,c4] + n1+n2;
% corners = floor(rand(4,1)*size(lrp.knots,1)+1+n1+n2)

if pressureType==4,
	e5 =    find(lrp.knots(: ,lrp.p(1)+1)<=0);
	e5 = e5(find(lrp.knots(e5,  end-1   )==0));
	e5 = e5(find(lrp.knots(e5,  end     )>0));
	e6 =    find(lrp.knots(: ,lrp.p(1)+1)==0);
	e6 = e6(find(lrp.knots(e6,  end-1   )<=0));
	e6 = e6(find(lrp.knots(e6,lrp.p(1)+2)>0 ));

	c1 = intersect(e1,e5);
	c5 = intersect(e5,e6);
	c6 = intersect(e6,e3);
	corners = [c1,c2,c3,c4,c5,c6] + n1+n2
	A(corners,:) = 0 ;
	A(:,corners) = 0 ;
	A(corners,corners) = eye(numel(corners));
	b(corners) = 0;
	b(c3) = -10;
	b(c4) = -10;
	return;
end

if pressureType>2
% 	if exist('pex')==1
% 		c1 = c1+n1+n2;
% 		c2 = c2+n1+n2;
% 		c3 = c3+n1+n2;
% 		c4 = c4+n1+n2;
% 		b       = b - A(:,c1)*pex(0,0);
% 		b(c1)   = pex(0,0);
% 
% 		b       = b - A(:,c2)*pex(0,1);
% 		b(c2)   = pex(0,1);
% 
% 		b       = b - A(:,c3)*pex(1,0);
% 		b(c3)   = pex(1,0);
% 
% 		b       = b - A(:,c4)*pex(1,1);
% 		b(c4)   = pex(1,1);
% 	else
  	b(corners)   = 0;
% 	end
	A(corners,:) = 0 ;
	A(:,corners) = 0 ;
	A(corners,corners) = eye(numel(corners));
end

if pressureType>1
	A = [A; zeros(1,n1+n2), avg_p']; % add the average zero preassure
	b = [b; 0];
end

