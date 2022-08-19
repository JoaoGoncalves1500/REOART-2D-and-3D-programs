function [vxy dx_vxy dy_vxy dx_dx_vxy dy_dx_vxy dx_dy_vxy...
    dy_dy_vxy] = cv2dr(v);
%===================================================================%
% This function interpolates the sound speed and the spatial 
% gradients with the position of the acoustic ray. This allows to 
% determinate the sound speed along the acoustic ray trajectory.
%
% The global variables are the sound speed and spatial gradients
%  matrices, distance and depth of the 2D simulation plan.
%===================================================================%
global Mvxy Mdx_vxy Mdy_vxy Mdx_dx_vxy Mdy_dx_vxy Mdy_dy_vxy...
    Mdx_dy_vxy dis_m depth

x=v(1); y=v(2); z=v(3); w=v(4);

[X,Y]=meshgrid(dis_m,depth);

vxy =interp2(X,Y,Mvxy,x,y,'linear',1500);

dx_vxy=interp2(X,Y,Mdx_vxy,x,y,'linear',0.005);

dy_vxy=interp2(X,Y,Mdy_vxy,x,y,'linear',0.0005);

dx_dx_vxy=interp2(X,Y,Mdx_dx_vxy,x,y,'linear',0.0005);

dy_dx_vxy=interp2(X,Y,Mdy_dx_vxy,x,y,'linear',0.0005);

dy_dy_vxy=interp2(X,Y,Mdy_dy_vxy,x,y,'linear',0.0005);

dx_dy_vxy=interp2(X,Y,Mdx_dy_vxy,x,y,'linear',0.0005);
end
    




