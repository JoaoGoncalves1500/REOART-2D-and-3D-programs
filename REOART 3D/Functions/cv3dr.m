function [vabc da_vabc db_vabc dc_vabc da_da_vabc db_da_vabc...
    dc_da_vabc da_db_vabc db_db_vabc dc_db_vabc da_dc_vabc...
    db_dc_vabc dc_dc_vabc dt_vabc da_dt_vabc db_dt_vabc...
    dc_dt_vabc] = cv3dr(v);
%===================================================================%
% This function interpolates the sound speed and the spatial 
% gradients with the position of the acoustic ray. This allows to 
% determinate the sound speed along the acoustic ray trajectory.
%
% The global variables are the sound speed and spatial 
% gradients matrices, distance and depth of the 3D simulation plan.
%===================================================================%
global Mvxyz Mdx_vxyz Mdy_vxyz Mdz_vxyz Mdx_dx_vxyz Mdx_dy_vxyz...
    Mdx_dz_vxyz Mdy_dx_vxyz Mdy_dy_vxyz Mdy_dz_vxyz Mdz_dx_vxyz...
    Mdz_dy_vxyz Mdz_dz_vxyz X_2 Y_2 Lat_3D Long_3D X_m Y_m depth...
    Depth_3D

a=v(1); b=v(2); c=v(3); x=v(4); y=v(5); z=v(6);

[X,Y,Z] = meshgrid(X_m(:,1)',Y_m(1,:),depth);

vabc=interp3(X,Y,Z,Mvxyz,a,b,c,'linear',1495);

da_vabc=interp3(X,Y,Z,Mdx_vxyz,a,b,c,'linear',0.0005);

db_vabc=interp3(X,Y,Z,Mdy_vxyz,a,b,c,'linear',0.0005);

dc_vabc=interp3(X,Y,Z,Mdz_vxyz,a,b,c,'linear',0.0005);

da_da_vabc=interp3(X,Y,Z,Mdx_dx_vxyz,a,b,c,'linear',0.0005);

db_da_vabc=interp3(X,Y,Z,Mdy_dx_vxyz,a,b,c,'linear',0.0005);

dc_da_vabc=interp3(X,Y,Z,Mdz_dx_vxyz,a,b,c,'linear',0.0005);

da_db_vabc=interp3(X,Y,Z,Mdx_dy_vxyz,a,b,c,'linear',0.0005);

db_db_vabc=interp3(X,Y,Z,Mdy_dy_vxyz,a,b,c,'linear',0.0005);

dc_db_vabc=interp3(X,Y,Z,Mdz_dy_vxyz,a,b,c,'linear',0.0005);

da_dc_vabc=interp3(X,Y,Z,Mdx_dz_vxyz,a,b,c,'linear',0.0005);

db_dc_vabc=interp3(X,Y,Z,Mdy_dz_vxyz,a,b,c,'linear',0.0005);

dc_dc_vabc=interp3(X,Y,Z,Mdz_dz_vxyz,a,b,c,'linear',0.0005);


dt_vabc=da_vabc*x+db_vabc*y+dc_vabc*z;
%
da_dt_vabc=da_da_vabc*x+da_db_vabc*y+da_dc_vabc*z;
db_dt_vabc=db_da_vabc*x+db_db_vabc*y+db_dc_vabc*z;
dc_dt_vabc=dc_da_vabc*x+dc_db_vabc*y+dc_dc_vabc*z;
end

