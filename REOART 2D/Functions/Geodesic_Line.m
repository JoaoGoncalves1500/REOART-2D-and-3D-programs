function [y_lat,x_long,dis_m,prof_area,long_area,lat_area,...
    inc_long,inc_lat]= Geodesic_Line(latitude,longitude,...
    elevation,lat_S,lat_N,long_W,long_E,lat_pos,long_pos,...
    azimute, distancia,inc_dis)
%===================================================================%
% This function calculates the geodesic line.
% The output are the geographical positions that compose the geodesic
% line.
%===================================================================%


inc_long=30/3600;
inc_lat=inc_long;
prof_max=min(min(elevation));
prof_min=max(max(elevation));
ind_long=round((long_pos-min(longitude))/inc_long+1);
ind_lat=round((lat_pos-min(latitude))/inc_lat+1);
prof_pos=elevation(ind_long,ind_lat);
ind_long_W=round((long_W-min(longitude))/inc_long+1);
ind_long_E=round((long_E-min(longitude))/inc_long+1);
ind_lat_S=round((lat_S-min(latitude))/inc_lat+1);
ind_lat_N=round((lat_N-min(latitude))/inc_lat+1);
prof_area=elevation(ind_long_W:ind_long_E,ind_lat_S:ind_lat_N);
size(prof_area);
max_prof_area=min(min(prof_area));
min_prof_area=max(max(prof_area));
long_area=long_W:inc_long:long_E;
size(long_area);
lat_area=lat_S:inc_lat:lat_N;
size(lat_area);
%-------------------------------------------------------------------%
distancia_m=distancia*1852; 
inc_dis_m=inc_dis*1852; 
dis_m=0:inc_dis_m:distancia_m; 
long_pos=-long_pos; 
s=dis_m; 
exc=0.0819919; 
a_eli=6378388; 
if azimute<=180
    azimute_S=180+azimute;
else
    azimute_S=azimute-180;
end
R_mer=(a_eli*(1-exc^2))/((1-exc^2*(sind(lat_pos))^2)^(3/2));
N=a_eli/((1-exc^2*sind(lat_pos)^2)^(1/2));
B=1/(R_mer*sind(1/3600));
C=tand(lat_pos)/(2*R_mer*N*sind(1/3600));
D=(3*exc^2*sind(lat_pos)*cosd(lat_pos)*sind(1/3600))/(2*(1-exc^2*...
    (sind(lat_pos))^2));
E=(1+3*(tand(lat_pos))^2)/(6*N^2);
del_fi=-s.*cosd(azimute_S).*B-s.^2.*(sind(azimute_S)).^2.*C;
h=s.*cosd(azimute_S).*B;
delta_fi_2linha=-s.*cosd(azimute_S).*B-s.^2.*(sind(azimute_S)).^2.*...
    C-del_fi.^2.*D+h.*s.^2.*(sind(azimute_S)).^2.*E;
delta_fi=delta_fi_2linha./3600;
lat_2=lat_pos+delta_fi;
N_linha=a_eli./((1-exc.^2.*sind(lat_2).^2).^(1/2));
A_linha=1./(N_linha.*sind(1./3600));
log_delta_lambda_2linha=log10(s)+log10(abs(sind(azimute_S)))+...
    log10(abs(A_linha))+log10(abs(secd(lat_2)));
delta_long=10.^log_delta_lambda_2linha./3600;
if azimute_S <= 180
    long_2=long_pos+delta_long;
else
    long_2=long_pos-delta_long;
end
y_lat=lat_2; 
x_long=-long_2; 
end

