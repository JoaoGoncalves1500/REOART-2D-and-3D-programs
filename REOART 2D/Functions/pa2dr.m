function [f,k,m]=pa2dr(t,v,vd)
%===================================================================%
% t - time
% v = [x;y;z;w]
% xd - derivative of x
% yd - derivative of y
% zd - derivative of z
% wd - derivative of w
% We suppose: xd^2+yd^2=1
%===================================================================%
x=v(1); y=v(2); z=v(3); w=v(4); 
% 
xd=vd(1); yd=vd(2); zd=vd(3); wd=vd(4);

% 
[vxy dx_vxy dy_vxy dx_dx_vxy dy_dx_vxy dx_dy_vxy dy_dy_vxy] = cv2dr(v);
% vxy- sound velocity field

% 
k_31=(w*dx_dx_vxy-z*dx_dy_vxy)/vxy-(w*dx_vxy-z*dy_vxy)*dx_vxy/vxy^2;
k_32=(w*dy_dx_vxy-z*dy_dy_vxy)/vxy-(w*dx_vxy-z*dy_vxy)*dy_vxy/vxy^2;
k_33=-wd-dy_vxy/vxy;
k_34=zd+dx_vxy/vxy;

k_41=0;
k_42=0;
k_43=zd;
k_44=wd;

% System of differential equations without equal sign and all reduced 
% to the left hand side
f=[xd-z;   
   yd-w;
   zd*w-wd*z+(w*dx_vxy-z*dy_vxy)/vxy;
  % zd*z+wd*w-vxy*(dx_vxy*z+dy_vxy*w)];
   zd*z+wd*w];

% Jacobian with respect to vector v: x,y,z,w
k=[0,0,-1,0;
   0,0,0,-1;
    k_31,k_32,k_33,k_34;
    k_41,k_42,k_43,k_44];

%Jacobian with respect to vector vd: xd, yd, zd, wd
m=[1,0,0,0;
   0,1,0,0;
   0,0,w,-z;
   0,0,z,w];
end
