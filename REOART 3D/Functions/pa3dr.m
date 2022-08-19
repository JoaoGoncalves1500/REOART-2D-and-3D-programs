function [f,k,m]=pa3dr(t,v,vd)
%===================================================================%
% t - time
% v = [a;b;c;x;y;z]
% ad - derivative of a
% bd - derivative of b
% cd - derivative of c
% xd - derivative of x
% yd - derivative of y
% zd - derivative of z
%===================================================================%
a=v(1); b=v(2); c=v(3); x=v(4); y=v(5); z=v(6);
%
ad=vd(1); bd=vd(2); cd=vd(3); xd=vd(4); yd=vd(5); zd=vd(6);

%
[vabc da_vabc db_vabc dc_vabc da_da_vabc db_da_vabc dc_da_vabc da_db_vabc...
    db_db_vabc dc_db_vabc da_dc_vabc db_dc_vabc dc_dc_vabc ...
    dt_vabc da_dt_vabc db_dt_vabc dc_dt_vabc]=cv3dr(v);
% vabc- sound velocity field

% System of differential equations without equal sign and all reduced 
% to the left hand side
f=[ad-x;   
   bd-y;
   cd-z;
   (x*z-(y^2+z^2))*xd+(x*y+y*z)*yd+(x*z-(x^2+y^2))*zd+(dt_vabc*(x+z)-...
   da_vabc-dc_vabc)/vabc;
   (x*y)*xd-(x^2+z^2)*yd+(y*z)*zd+(dt_vabc*y-db_vabc)/vabc;
   x*xd+y*yd+z*zd];

k_41=(da_dt_vabc*(x+z)-da_da_vabc-da_dc_vabc)/vabc-da_vabc*...
    (dt_vabc*(x+z)-da_vabc-dc_vabc)/vabc^2;
%
k_42=(db_dt_vabc*(x+z)-db_da_vabc-db_dc_vabc)/vabc-db_vabc*...
    (dt_vabc*(x+z)-da_vabc-dc_vabc)/vabc^2;
%
k_43=(dc_dt_vabc*(x+z)-dc_da_vabc-dc_dc_vabc)/vabc-dc_vabc*...
    (dt_vabc*(x+z)-da_vabc-dc_vabc)/vabc^2;
%
k_44=z*xd+y*yd+z*zd+dt_vabc/vabc;
%
k_45=-2*y*xd+x*yd-2*y*zd;
%
k_46=x-2*z*xd+y*yd+x*zd+dt_vabc*x/vabc;
%
%
k_51=(da_dt_vabc*y-da_db_vabc)/vabc-da_vabc*(dt_vabc*y-db_vabc)/vabc^2;
k_52=(db_dt_vabc*y-db_db_vabc)/vabc-db_vabc*(dt_vabc*y-db_vabc)/vabc^2;
k_53=(dc_dt_vabc*y-dc_db_vabc)/vabc-dc_vabc*(dt_vabc*y-db_vabc)/vabc^2;
%
k_54=y*xd-2*x*yd;
k_55=x*xd+z*zd+dt_vabc/vabc;
k_56=-2*z*yd+y*zd;
%
k_61=0;
k_62=0;
k_63=0;
%
k_64=xd;
k_65=yd;
k_66=zd;

% Jacobian with respect to vector v: a,b,c,x,y,z
k=[0,0,0,-1,0,0;
   0,0,0,0,-1,0;
   0,0,0,0,0,-1;
    k_41,k_42,k_43,k_44,k_45,k_46;
    k_51,k_52,k_53,k_54,k_55,k_56;
    k_61,k_62,k_63,k_64,k_65,k_66];

%Jacobian with respect to vector vd: ad, bd, cd, xd, yd, zd
m=[1,0,0,0,0,0;
   0,1,0,0,0,0;
   0,0,1,0,0,0;
   0,0,0,(x*z-(y^2+z^2)),(x*y+y*z),(x*z-(x^2+y^2));
   0,0,0,x*y,-(x^2+z^2),y*z;
   0,0,0,x,y,z];
end




