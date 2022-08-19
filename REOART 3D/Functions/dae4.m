function ys=dae4(f,tspan,y0,nint,g)
% function ys=dae4(f,tspan,y0,nint,g)
% solves a set of differential algebraic equations (DAEs)
%        f(t,y,y')=0  where y'=dy/dt
% with a 4th order method starting from y0 at time t0 and 
% finishing at time tfin where tspan=[t0 t1 ... tfin].
% y0 is a column vector, tspan is a row vector.
%
% The solution is returned at all the times in tspan.
% The time steps are approx diff(ts)/nint within the domain.
% Error management is entirely up to the user via tspan & nint.
% If optional nint is omitted, then it is assumed to be 1.
% If matlab warms of matrices with high condition number,
% then increase nint.
%
% The jacobians of f, namely k=df/dy and m=df/dy', must be 
% provided by f, at each time compute: [f,k,m]=func(t,y,y').
% Both m and k may be returned as sparse matrices.
% If warnings of poor convergence occur, then the coded
% jacobians probably have errors.
%
% The optional argument g is the name of a user supplied
% function g(t,y,y') that is invoked immediately the
% solution is computed at the times in tspan.
%
% The initial state y0 should be consistent with the 
% algebraic part of the DAE, but if not consistent then 
% transient oscillations will appear as it works towards
% consistency.
%
% The method will also work well for stiff sets of ODEs.
% It is unsuitable for problems with an oscillatory spectrum
% as non-dissipative oscillations will grow for h*omega<4 ; 
% instead use dae4o.m
%
% See pendrun.m, penddae.m & pendg.m for a pendulum example.
% See also dae2.m and dae4o.m for other versions.
%
% (c) Tony Roberts, 18 Aug 1998, aroberts@usq.edu.au

if nargin<4, nint=1; end
gcall=(nargin==5);
nout=length(tspan);
ndim=length(y0);
ys=zeros(ndim,nout);
newtol=1e-6;
newtmax=10;
newtit=zeros(1,newtmax);
% weights for BDF est of deriv
wd=[0 1 1/2 1/3 1/4]; 
wds=sum(wd);

% to allow for varying spaced output, fit a spline
% and solve in s=[1,nts] rather than in t
% dt is dt/ds at each time s
nts=1+nint*(nout-1);
if nts<5, disp('ERROR: dae4 needs at least 4 steps'), return, end
ts=spline(1:nint:nts,tspan,(1:nts));
dt=spline(1:nint:nts,tspan,(1:nts)+1e-7);
dt=(dt-ts)/1e-7;

% initialise by solving first four steps together assuming
% a quartic between them (zero fifth difference).
yy=zeros(ndim,4); % initial guess
% y0
% yy
y=[y0 yy]; 
for newt=1:newtmax
	% extrapolate
	y2=rot90(cumsum(rot90(y )),-1);
	y3=rot90(cumsum(rot90(y2)),-1);
	y4=rot90(cumsum(rot90(y3)),-1);
	y5=rot90(cumsum(rot90(y4)),-1);
	% evaluate residuals
	[f2,k2,m2]=feval(f,ts(2),y2(:,1),y2*wd'/dt(2));
	[f3,k3,m3]=feval(f,ts(3),y3(:,1),y3*wd'/dt(3));
	[f4,k4,m4]=feval(f,ts(4),y4(:,1),y4*wd'/dt(4));
	[f5,k5,m5]=feval(f,ts(5),y5(:,1),y5*wd'/dt(5));
	% solve simultaneous equations
	yy(:)=-[  k2+m2/dt(2)    k2+3/2*m2/dt(2)    k2+11/6*m2/dt(2)    k2+25/12*m2/dt(2)
	        2*k3+m3/dt(3)  3*k3+5/2*m3/dt(3)  4*k3+13/3*m3/dt(3)  5*k3+77/12*m3/dt(3)
	        3*k4+m4/dt(4)  6*k4+7/2*m4/dt(4) 10*k4+47/6*m4/dt(4) 15*k4+171/12*m4/dt(4)
	        4*k5+m5/dt(5) 10*k5+9/2*m5/dt(5) 20*k5+37/3*m5/dt(5) 35*k5+319/12*m5/dt(5)
			]\[f2;f3;f4;f5];
	y(:,2:5)=y(:,2:5)+yy;
	if max(abs(yy(:)))<newtol*max(abs(y(:))), break, end
end
% output as requested
ys(:,1)=y(:,1);
if gcall, feval(g,ts(1),y(:,1),y*wd'/dt(1)); end
for n=2:5
	y=rot90(cumsum(rot90(y)),-1);
    if rem(n-1,nint)==0, ys(:,1+(n-1)/nint)=y(:,1); 
    if gcall, feval(g,ts(n),y(:,1),y*wd'/dt(n)); end, end
end

% take fourth order steps over domain
for n=6:nts
	h=dt(n);
	y=rot90(cumsum(rot90(y)),-1); % extrapolate a guess
	for newt=1:newtmax
    	[f1,k1,m1]=feval(f,ts(n),y(:,1),y*wd'/h);
    	w=-(wds*m1/h+k1)\f1;
    	y=y+w*ones(size(wd));
    	if max(abs(w))<newtol*max(abs(y(:))), break, end
	end
	newtit(newt)=newtit(newt)+1;
    if rem(n-1,nint)==0, ys(:,1+(n-1)/nint)=y(:,1); 
    if gcall, feval(g,ts(n),y(:,1),y*wd'/dt(n)); end, end
end

% check on how many Newtonian iterations were required
newtm=sum((1:newtmax).*newtit)/sum(newtit);
if newtm>newtmax/2,
   disp('WARNING: poor or no convergence in dae4')
   newtit=newtit
end
