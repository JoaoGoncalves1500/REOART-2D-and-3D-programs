%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------%
%                             REOART 3D 
%-------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REOART 2D (Real Environment Ocean-Acoustic Ray Tracing) is an 
% ocean-acoustic ray tracing model with the ability to study sound 
% propagation in a real ocean environment. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this script is where the sound speed and bathymetry are 
% calculated 
%-------------------------------------------------------------------%
clearvars -global; close all; clear; clc;
global Mvxyz Mdx_vxyz Mdy_vxyz Mdz_vxyz Mdx_dx_vxyz Mdx_dy_vxyz...
   Mdx_dz_vxyz Mdy_dx_vxyz Mdy_dy_vxyz Mdy_dz_vxyz Mdz_dx_vxyz...
   Mdz_dy_vxyz Mdz_dz_vxyz X_2 Y_2 Lat_3D Long_3D X_m Y_m...
   depth Depth_3D prof_area

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------Initial geographic position---------------------%
lat_pos=41;
long_pos=-8.8;
azimute = 270;
azimute_1 = azimute+90;
distance = 50.8;
inc_dis = 0.5; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------%
%                        GEBCO DATABASE                             % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Selecting the area of the world. 
% In this case is just the North Atlantic
% area (gebco=1).
gebco = 1;  
if gebco==1
    longitude=ncread('GEBCO_AtlanticoNorte.nc','lon');
    latitude=ncread('GEBCO_AtlanticoNorte.nc','lat');
    elevation=ncread('GEBCO_AtlanticoNorte.nc','elevation');
end

% Central position geodesic
[y_lat,x_long,dis_m]=Geodesica_GEBCO(lat_pos,long_pos,azimute,...
    distance,inc_dis); 

% position 1 geodesic
[y_lat_1,x_long_1,dis_m_1]=Geodesica_GEBCO(lat_pos,long_pos,azimute_1,...
    distance,inc_dis); 

% Central position 
posicao_geo = [x_long;y_lat];
posicao_geo=posicao_geo';

% Position 1 
posicao_geo_1 = [x_long_1;y_lat_1];
posicao_geo_1=posicao_geo_1';

[X]=meshgrid(posicao_geo(:,1),posicao_geo_1(:,1) ); % Longitude
[Y]=meshgrid(posicao_geo_1(:,2),posicao_geo(:,2)); % Latitude
X=X';

% Positions in meters
[X_m]=meshgrid(dis_m(1,:),dis_m(1,:));
[Y_m]=meshgrid(dis_m_1(1,:),dis_m_1(1,:));
X_m=X_m';

X_2=X(:,1); %Mx1 matrix - Longitude
Y_2=Y(1,:); %Mx1 matrix - Latitude

area=8; % Change the number to change the simulation area 
if area==1
    lat_S=34.5; lat_N=44.5; long_W=-15.5; long_E=-7.5; % North 
elseif area==2                                         % Atlantic
    lat_S=30; lat_N=45; long_W=-40; long_E=-5; % Northeast 
elseif area==3                                 % Atlantic
    lat_S=28.5; lat_N=36.5; long_W=-20.5; long_E=-13.5; % Madeira
elseif area==4
    lat_S=36; lat_N=40; long_W=-33; long_E=-22; % Açores
elseif area==5
    lat_S=34; lat_N=44; long_W=-15; long_E=-6; % Continental Portugal
elseif area==6
     lat_S=30; lat_N=43; long_W=-32; long_E=-6; % Portugal ZEE
elseif area==7
    lat_S=46; lat_N=61; long_W=-25; long_E=0; % Scotland
elseif area==8
    lat_S=min(min(Y)); lat_N=max(max(Y)); long_W=min(min(X)); long_E=max(max(X));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------%
%                        GEODESIC LINE                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[y_lat,x_long,dis_m,prof_area,long_area,lat_area,inc_long,...
    inc_lat]= Geodesic_Line(latitude,longitude,elevation,lat_S,lat_N,...
    long_W,long_E,lat_pos,long_pos,azimute, distance,inc_dis);

ind_x_long=round((x_long-min(longitude))./inc_long+1);
ind_y_lat=round((y_lat-min(latitude))./inc_lat+1);
perfil_fundo(1:length(ind_x_long))=0;
for i=1:length(ind_x_long)
perfil_fundo(i)=elevation(ind_x_long(i),ind_y_lat(i));
end
perfil_fundo;
dis_km=dis_m./1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------%
%                        TEOS-10 FUNCTIONS                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----- WOD18 Files
%loading the salinity and temperature values from WOD18 database
WOD18_T=csvread('woa18_95A4_t00an04.csv',2,0);
WOD18_S=csvread('woa18_95A4_s00an04.csv',2,0);

linhas_WOD18_T=size(WOD18_T,1);
colunas_WOD18_T=size(WOD18_T,2); 
linhas_WOD18_S=size(WOD18_S,1);
colunas_WOD18_S=size(WOD18_S,2);

lat_T=WOD18_T(1:linhas_WOD18_T,1);
long_T=WOD18_T(1:linhas_WOD18_T,2);
T(1:linhas_WOD18_T,1:colunas_WOD18_T-2)=WOD18_T(1:linhas_WOD18_T,3:...
    colunas_WOD18_T);
lat_S=WOD18_S(1:linhas_WOD18_S,1);
long_S=WOD18_S(1:linhas_WOD18_S,2);
S(1:linhas_WOD18_S,1:colunas_WOD18_S-2)=WOD18_S(1:linhas_WOD18_S,3:...
    colunas_WOD18_S);

%Salinity and temperature matrices
if area==7 % For Scotland simulation area
    % Salinity
    posicao_geoS_1=[lat_S long_S];
    posicao_geoS=posicao_geoS_1(478830:585857,1:2);
    S=S(478830:585857,:);
    % Temperature
    posicao_geoT_1=[lat_T long_T];
    posicao_geoT=posicao_geoS;
    T=T(478830:585857,:);
 else % All other simulation area
    % Salinity
    posicao_geoS_1=[lat_S long_S];
    posicao_geoS=posicao_geoS_1(478830:533094,1:2);
    S=S(478830:533094,:);
    % Temperature
    posicao_geoT_1=[lat_T long_T];
    posicao_geoT=posicao_geoS;
    T=T(478830:533094,:);
end

X_1=reshape(X,1,[]); %transform matrix X into a matrix with 1 column
Y_1=reshape(Y,1,[]); %transform matrix Y into a matrix with 1 column

posicao_geo3D = [Y_1;X_1];
posicao_geo3D = posicao_geo3D';

% Finding the nearest neighbors between data and geographic positions 
Idx_S = knnsearch(posicao_geoS(:,1:2),posicao_geo3D(:,1:2));
Idx_T = knnsearch(posicao_geoT(:,1:2),posicao_geo3D(:,1:2));

Idx_S = reshape(Idx_S,length(X),[]); 
Idx_T = reshape(Idx_T,length(X),[]);

% Depth values in meters
depth=[0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,...
    95,100,125,150,175,200,225,250,275,300,325,350,375,400,425,...
    450,475,500,550,600,650,700,750,800,850,900,950,1000,1050,...
    1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1600,1650,...
    1700,1750,1800,1850,1900,1950,2000,2100,2200,2300,2400,2500,...
    2600,2700,2800,2900,3000,3100,3200,3300,3400,3500,3600,3700,...
    3800,3900,4000,4100,4200,4300,4400,4500,4600,4700,4800,4900,...
    5000,5100,5200,5300,5400,5500];

% Construction of Salinity and Temperature matrices
for i=1:length(X)
    for ii=1:length(depth)
        for j=1:length(Y)
            Practical_Salinity(i,j,ii)=S(Idx_S(i),ii);
            Insitu_Temperature(i,j,ii)=T(Idx_T(i),ii);
        end 
    end
end

% Pressure calculation
for i=1:length(Y)
    for j=1:length(Y)
        p(i,j,1:length(depth)) = gsw_p_from_z(-depth,Y(i,j));
    end
end

% 3D Longitude
Long_3D=zeros(length(X),length(Y),length(depth));
for k=1:length(depth)
    Long_3D(:,:,k)=depth(1,k);
end
for j=1:length(Y)
    Long_3D(:,j,:)=Y(1,j);
end
for i=1:length(X)
    Long_3D(i,:,:)=X(i,1); 
end

% 3D Latitude
Lat_3D = zeros(length(X),length(Y),length(depth));
for k=1:length(depth)
    Lat_3D(:,:,k)=depth(1,k);
end
for i=1:length(X)
    Lat_3D(i,:,:)=X(i,1);
end
for j=1:length(Y)
    Lat_3D(:,j,:)=Y(1,j); 
end

%  3D Depth
Depth_3D = zeros(length(X),length(Y),length(depth));
for i=1:length(X)
    Depth_3D(i,:,:)=X(i,1);
end
for j=1:length(Y)
    Depth_3D(:,j,:)=Y(1,j);
end
for k=1:length(depth)
    Depth_3D(:,:,k)=depth(1,k); 
end

% 3D absolute salinity calculation
SA = gsw_SA_from_SP(Practical_Salinity,p,Long_3D,Lat_3D);    

% 3D conservative temperature calculation
CT = gsw_CT_from_t(SA,Insitu_Temperature,p);
 
%  3D sound speed calculation
sound_speed = gsw_sound_speed(SA,CT,p);
Mvxyz=sound_speed;

%-------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Calculating spatial gradients of the sound speed           %
%-------------------------------------------------------------------%
% dx_vxyz dy_vxyz dz_vxyz - First derivatives
% dx_dx_vxyz dx_dy_vxyz dx_dz_vxyz -   -|
% dy_dx_vxyz dy_dy_vxyz dy_dz_vxyz -    |   Second derivatives
% dz_dx_vxyz dz_dy_vxyz dz_dz_vxyz -   -|
%-------------------------------------------------------------------%
% x - Longitude
% y - Latitude
% z - Depth

%-----Centered finite-difference formulas----
% First derivatives: f'(xi)= -f(xi+2)+8f(xi+1)-8f(xi-1)+f(xi-2)/12h
% Second derivatives: f''(xi)= -f(xi+2)+16f(xi+1)-30f(xi)+...
%                               16f(xi-1)-f(xi-2)/12h^2
% h is the increment used between geographic positions

%-----Backward finite-difference formulas----
% First derivatives: f'(xi)= 3f(xi)-4f(xi-1)+f(xi-2)/2h
% Second derivatives: f''(xi)= 2f(xi)-5f(xi-1)+4f(xi-2)-f(xi-3)/h^2

%-----Forward finite-difference formulas-----
% First derivatives: f'(xi)= -f(xi+2)+4f(xi+1)-3f(xi)/2h
% Second derivatives: f''(xi)= -f(xi+3)+4f(xi+2)-5f(xi+1)+2f(xi)/h^2

%-------------------------------------------------------------------%
%%%%%%%%%%$%%%%%%%%%%%% First Derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%
% dx - variation along the lines (Depth)
% dy - variation along the columns (Latitude)
% dz - variation between matrices (Longitude)

imersoes_2=[0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,...
    95,100,125,150,175,200,225,250,275,300,325,350,375,400,425,...
    450,475,500,550,600,650,700,750,800,850,900,950,1000,1050,...
    1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1600,1650,...
    1700,1750,1800,1850,1900,1950,2000,2100,2200,2300,2400,2500,...
    2600,2700,2800,2900,3000,3100,3200,3300,3400,3500,3600,3700,...
    3800,3900,4000,4100,4200,4300,4400,4500,4600,4700,4800,4900,...
    5000,5100,5200,5300,5400,5500,5600]; 
imersoes_2=imersoes_2';

%--------- dx_vxyz calculation -----------
hx=inc_dis*1852;
for i = 1:length(depth);
    for k = 1:length(Y);
        for ii = 1:length(X);
            if i==1 % Forward finite-difference method
                Mdx_vxyz(i,ii,k)= (-Mvxyz(i+2,ii,k)+4.*...
                    Mvxyz(i+1,ii,k)-3.*Mvxyz(i,ii,k))/(2*hx); 
            elseif i==2
                Mdx_vxyz(i,ii,k)= (-Mvxyz(i+2,ii,k)+4.*...
                    Mvxyz(i+1,ii,k)-3.*Mvxyz(i,ii,k))/(2*hx);
            elseif i==length(depth)% Backward finite-difference method
                Mdx_vxyz(i,ii,k)= (3.*Mvxyz(i,ii,k)-4.*...
                    Mvxyz(i-1,ii,k)+Mvxyz(i-2,ii,k))/(2*hx);
            elseif i==(length(depth)-1)
                Mdx_vxyz(i,ii,k)= (3.*Mvxyz(i,ii,k)-4.*...
                    Mvxyz(i-1,ii,k)+Mvxyz(i-2,ii,k))/(2*hx);
            else % Centered finite-difference method
                Mdx_vxyz(i,ii,k)= (-Mvxyz(i+2,ii,k)+8.*...
                    Mvxyz(i+1,ii,k)-8.*Mvxyz(i-1,ii,k)+...
                    Mvxyz(i-2,ii,k))/(12*hx);                                
            end
        end
    end
   
end

%--------- dy_vxyz calculation ---------------
hy=inc_dis*1852;
for i = 1:length(depth);
    for k = 1:length(Y);
        for ii = 1:length(X);
            if ii==1 %  Forward finite-difference method
                Mdy_vxyz(i,ii,k)= (-Mvxyz(i,ii+2,k)+4.*...
                    Mvxyz(i,ii+1,k)-3.*Mvxyz(i,ii,k))/(2*hy); 
            elseif ii==2
                Mdy_vxyz(i,ii,k)= (-Mvxyz(i,ii+2,k)+4.*...
                    Mvxyz(i,ii+1,k)-3.*Mvxyz(i,ii,k))/(2*hy); 
            elseif ii==length(X)% Backward finite-difference method
                Mdy_vxyz(i,ii,k)= (3.*Mvxyz(i,ii,k)-4.*...
                    Mvxyz(i,ii-1,k)+Mvxyz(i,ii-2,k))/(2*hy);
            elseif ii==(length(X)-1)
                Mdy_vxyz(i,ii,k)= (3.*Mvxyz(i,ii,k)-4.*...
                    Mvxyz(i,ii-1,k)+Mvxyz(i,ii-2,k))/(2*hy);
            else % Centered finite-difference method
                Mdy_vxyz(i,ii,k)= (-Mvxyz(i,ii+2,k)+8.*...
                    Mvxyz(i,ii+1,k)-8.*Mvxyz(i,ii-1,k)+...
                    Mvxyz(i,ii-2,k))/(12*hy);                                
            end
        end
    end
end

%--------- dz_vxyz calculation --------------
hz=diff(imersoes_2);
for k =1:length(depth)
   for i = 1:length(X);
       for ii = 1:length(Y);
           if ii==1 %  Forward finite-difference method
               Mdz_vxyz(k,i,ii)= (-Mvxyz(k,i,ii+2)+4.*...
                   Mvxyz(k,i,ii+1)-3.*Mvxyz(k,i,ii))/(2*hz(i,1)); 
           elseif ii==2
               Mdz_vxyz(k,i,ii)= (-Mvxyz(k,i,ii+2)+4.*...
                   Mvxyz(k,i,ii+1)-3.*Mvxyz(k,i,ii))/(2*hz(i,1)); 
           elseif ii==length(X)% Backward finite-difference method
               Mdz_vxyz(k,i,ii)= (3.*Mvxyz(k,i,ii)-4.*...
                   Mvxyz(k,i,ii-1)+Mvxyz(k,i,ii-2))/(2*hz(i,1));
           elseif ii==(length(X)-1)
               Mdz_vxyz(k,i,ii)= (3.*Mvxyz(k,i,ii)-4.*...
                   Mvxyz(k,i,ii-1)+Mvxyz(k,i,ii-2))/(2*hz(i,1));
           else % Centered finite-difference method
               Mdz_vxyz(k,i,ii)= (-Mvxyz(k,i,ii+2)+8.*...
                   Mvxyz(k,i,ii+1)-8.*Mvxyz(k,i,ii-1)+...
                   Mvxyz(k,i,ii-2))/(12*hz(i,1));                                
           end
       end
   end
end 

%-------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%% Second Derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%
% The derivatives: dx_dx_vxyz dx_dy_vxyz dx_dz_vxyz 
%                  dy_dx_vxyz dy_dy_vxyz dy_dz_vxyz     
%                  dz_dx_vxyz dz_dy_vxyz dz_dz_vxyz   
%-------------------------------------------------------------------%

%--------- dy_dx_vxyz calculation ---------
hy=inc_dis*1852;
for i = 1:length(depth);
    for k = 1:length(Y);
        for ii = 1:length(X);
            if ii==1 %  Forward finite-difference method
                Mdy_dx_vxyz(i,ii,k)= (-Mdx_vxyz(i,ii+2,k)+4.*...
                    Mdx_vxyz(i,ii+1,k)-3.*Mdx_vxyz(i,ii,k))/(2*hy); 
            elseif ii==2
                Mdy_dx_vxyz(i,ii,k)= (-Mdx_vxyz(i,ii+2,k)+4.*...
                    Mdx_vxyz(i,ii+1,k)-3.*Mdx_vxyz(i,ii,k))/(2*hy); 
            elseif ii==length(X)% Backward finite-difference method
                Mdy_dx_vxyz(i,ii,k)= (3.*Mdx_vxyz(i,ii,k)-4.*...
                    Mdx_vxyz(i,ii-1,k)+Mdx_vxyz(i,ii-2,k))/(2*hy);
            elseif ii==(length(X)-1)
                Mdy_dx_vxyz(i,ii,k)= (3.*Mdx_vxyz(i,ii,k)-4.*...
                    Mdx_vxyz(i,ii-1,k)+Mdx_vxyz(i,ii-2,k))/(2*hy);
            else % Centered finite-difference method
                Mdy_dx_vxyz(i,ii,k)= (-Mdx_vxyz(i,ii+2,k)+8.*...
                    Mdx_vxyz(i,ii+1,k)-8.*Mdx_vxyz(i,ii-1,k)+...
                    Mdx_vxyz(i,ii-2,k))/(12*hy);                                
            end
        end 
    end
end

%--------- dx_dx_vxyz calculation --------- 
hx=inc_dis*1852;
for i = 1:length(depth);
    for k = 1:length(Y);
        for ii = 1:length(X);
            if i==1 %  Forward finite-difference method
                Mdx_dx_vxyz(i,ii,k)= (-Mdx_vxyz(i+2,ii,k)+4.*...
                    Mdx_vxyz(i+1,ii,k)-3.*Mdx_vxyz(i,ii,k))/(2*hx); 
            elseif i==2
                Mdx_dx_vxyz(i,ii,k)= (-Mdx_vxyz(i+2,ii,k)+4.*...
                    Mdx_vxyz(i+1,ii,k)-3.*Mdx_vxyz(i,ii,k))/(2*hx);
            elseif i==length(depth)% Backward finite-difference method
                Mdx_dx_vxyz(i,ii,k)= (3.*Mdx_vxyz(i,ii,k)-4.*...
                    Mdx_vxyz(i-1,ii,k)+Mdx_vxyz(i-2,ii,k))/(2*hx);
            elseif i==(length(depth)-1)
                Mdx_dx_vxyz(i,ii,k)= (3.*Mdx_vxyz(i,ii,k)-4.*...
                    Mdx_vxyz(i-1,ii,k)+Mdx_vxyz(i-2,ii,k))/(2*hx);
            else % Centered finite-difference method
                Mdx_dx_vxyz(i,ii,k)= (-Mdx_vxyz(i+2,ii,k)+8.*...
                    Mdx_vxyz(i+1,ii,k)-8.*Mdx_vxyz(i-1,ii,k)+...
                    Mdx_vxyz(i-2,ii,k))/(12*hx);                                
            end
        end
    end
end

%--------- dz_dx_vxyz calculation --------- 
hz=diff(imersoes_2);
for k = 1:length(depth)
    for i = 1:length(X);
        for ii = 1:length(Y);
            if ii==1 % Forward finite-difference method
                Mdz_dx_vxyz(k,i,ii)= (-Mdx_vxyz(k,i,ii+2)+4.*...
                    Mdx_vxyz(k,i,ii+1)-3.*Mdx_vxyz(k,i,ii))/...
                    (2*hz(i,1)); 
            elseif ii==2
                Mdz_dx_vxyz(k,i,ii)= (-Mdx_vxyz(k,i,ii+2)+4.*...
                    Mdx_vxyz(k,i,ii+1)-3.*Mdx_vxyz(k,i,ii))/...
                    (2*hz(i,1)); 
            elseif ii==length(X)% Backward finite-difference method
                Mdz_dx_vxyz(k,i,ii)= (3.*Mdx_vxyz(k,i,ii)-4.*...
                    Mdx_vxyz(k,i,ii-1)+Mdx_vxyz(k,i,ii-2))/...
                    (2*hz(i,1));
            elseif ii==(length(X)-1)
                Mdz_dx_vxyz(k,i,ii)= (3.*Mdx_vxyz(k,i,ii)-4.*...
                    Mdx_vxyz(k,i,ii-1)+Mdx_vxyz(k,i,ii-2))/...
                    (2*hz(i,1));
            else % Centered finite-difference method
                Mdz_dx_vxyz(k,i,ii)= (-Mdx_vxyz(k,i,ii+2)+8.*...
                    Mdx_vxyz(k,i,ii+1)-8.*Mdx_vxyz(k,i,ii-1)+...
                    Mdx_vxyz(k,i,ii-2))/(12*hz(i,1));                                
            end
        end
    end
end

%--------- dy_dy_vxy calculation ---------
hy=inc_dis*1852;
for i = 1:length(depth);
    for k = 1:length(Y);
        for ii = 1:length(X);
            if ii==1 %  Forward finite-difference method
                Mdy_dy_vxyz(i,ii,k)= (-Mdy_vxyz(i,ii+2,k)+4.*...
                    Mdy_vxyz(i,ii+1,k)-3.*Mdy_vxyz(i,ii,k))/(2*hy); 
            elseif ii==2
                Mdy_dy_vxyz(i,ii,k)= (-Mdy_vxyz(i,ii+2,k)+4.*...
                    Mdy_vxyz(i,ii+1,k)-3.*Mdy_vxyz(i,ii,k))/(2*hy);
            elseif ii==length(X)% Backward finite-difference method
                Mdy_dy_vxyz(i,ii,k)= (3.*Mdy_vxyz(i,ii,k)-4.*...
                    Mdy_vxyz(i,ii-1,k)+Mdy_vxyz(i,ii-2,k))/(2*hy);
            elseif ii==(length(X)-1)
                Mdy_dy_vxyz(i,ii,k)= (3.*Mdy_vxyz(i,ii,k)-4.*...
                    Mdy_vxyz(i,ii-1,k)+Mdy_vxyz(i,ii-2,k))/(2*hy);
            else % Centered finite-difference method
                Mdy_dy_vxyz(i,ii,k)= (-Mdy_vxyz(i,ii+2,k)+8.*...
                    Mdy_vxyz(i,ii+1,k)-8.*Mdy_vxyz(i,ii-1,k)+...
                    Mdy_vxyz(i,ii-2,k))/(12*hy);                                
            end
        end
    end
end

%--------- dx_dy_vxy calculation ---------
hx=inc_dis*1852;
for i = 1:length(depth);
    for k = 1:length(Y);
        for ii = 1:length(X);
            if i==1 %  Forward finite-difference method
                Mdx_dy_vxyz(i,ii,k)= (-Mdy_vxyz(i+2,ii,k)+4.*...
                    Mdy_vxyz(i+1,ii,k)-3.*Mdy_vxyz(i,ii,k))/(2*hx); 
            elseif i==2
                Mdx_dy_vxyz(i,ii,k)= (-Mdy_vxyz(i+2,ii,k)+4.*...
                    Mdy_vxyz(i+1,ii,k)-3.*Mdy_vxyz(i,ii,k))/(2*hx);
            elseif i==length(depth)% Backward finite-difference method
                Mdx_dy_vxyz(i,ii,k)= (3.*Mdy_vxyz(i,ii,k)-4.*...
                    Mdy_vxyz(i-1,ii,k)+Mdy_vxyz(i-2,ii,k))/(2*hx);
            elseif i==(length(depth)-1)
                Mdx_dy_vxyz(i,ii,k)= (3.*Mdy_vxyz(i,ii,k)-4.*...
                    Mdy_vxyz(i-1,ii,k)+Mdy_vxyz(i-2,ii,k))/(2*hx);
            else % Centered finite-difference method
                Mdx_dy_vxyz(i,ii,k)= (-Mdy_vxyz(i+2,ii,k)+8.*...
                    Mdy_vxyz(i+1,ii,k)-8.*Mdy_vxyz(i-1,ii,k)+...
                    Mdy_vxyz(i-2,ii,k))/(12*hx);                                
            end
        end
    end
end

%--------- dz_dy_vxy calculation ---------
hz=diff(imersoes_2);
for i = 1:length(depth);
    for k = 1:length(X);
        for ii = 1:length(Y);
            if ii==1 %  Forward finite-difference method
                Mdz_dy_vxyz(i,k,ii)= (-Mdy_vxyz(i,k,ii+2)+4.*...
                    Mdy_vxyz(i,k,ii+1)-3.*Mdy_vxyz(i,k,ii))/...
                    (2*hz(i,1)); 
            elseif ii==2
                Mdz_dy_vxyz(i,k,ii)= (-Mdy_vxyz(i,k,ii+2)+4.*...
                    Mdy_vxyz(i,k,ii+1)-3.*Mdy_vxyz(i,k,ii))/...
                    (2*hz(i,1)); 
            elseif ii==length(Y)% Backward finite-difference method
                Mdz_dy_vxyz(i,k,ii)= (3.*Mdy_vxyz(i,k,ii)-4.*...
                    Mdy_vxyz(i,k,ii-1)+Mdy_vxyz(i,k,ii-2))/...
                    (2*hz(i,1));
            elseif ii==(length(Y)-1)
                Mdz_dy_vxyz(i,k,ii)= (3.*Mdy_vxyz(i,k,ii)-4.*...
                    Mdy_vxyz(i,k,ii-1)+Mdy_vxyz(i,k,ii-2))/...
                    (2*hz(i,1));
            else % Centered finite-difference method
                Mdz_dy_vxyz(i,k,ii)= (-Mdy_vxyz(i,k,ii+2)+8.*...
                    Mdy_vxyz(i,k,ii+1)-8.*Mdy_vxyz(i,k,ii-1)+...
                    Mdy_vxyz(i,k,ii-2))/(12*hz(i,1));                                
            end
        end
    end
end

%--------- dx_dz_vxy calculation ---------
hx=inc_dis*1852;
for k = 1:length(depth)
    for i = 1:length(X);
        for ii = 1:length(Y);
            if k==1 %  Forward finite-difference method
                Mdx_dz_vxyz(k,i,ii)= (-Mdz_vxyz(k+2,i,ii)+4.*...
                    Mdz_vxyz(k+1,i,ii)-3.*Mdz_vxyz(k,i,ii))/(2*hx); 
            elseif k==2
                Mdx_dz_vxyz(k,i,ii)= (-Mdz_vxyz(k+2,i,ii)+4.*...
                    Mdz_vxyz(k+1,i,ii)-3.*Mdz_vxyz(k,i,ii))/(2*hx); 
            elseif k==length(depth)% Backward finite-difference method
                Mdx_dz_vxyz(k,i,ii)= (3.*Mdz_vxyz(k,i,ii)-4.*...
                    Mdz_vxyz(k-1,i,ii)+Mdz_vxyz(k-2,i,ii))/(2*hx);
            elseif k==(length(depth)-1)
                Mdx_dz_vxyz(k,i,ii)= (3.*Mdz_vxyz(k,i,ii)-4.*...
                    Mdz_vxyz(k-1,i,ii)+Mdz_vxyz(k-2,i,ii))/(2*hx);
            else % Centered finite-difference method
                Mdx_dz_vxyz(k,i,ii)= (-Mdz_vxyz(k+2,i,ii)+8.*...
                    Mdz_vxyz(k+1,i,ii)-8.*Mdz_vxyz(k-1,i,ii)+...
                    Mdz_vxyz(k-2,i,ii))/(12*hx);                                
            end
        end
    end
end

%--------- dy_dz_vxy calculation ---------
hy=inc_dis*1852;
for i = 1:length(depth);
    for k = 1:length(Y);
        for ii = 1:length(X);
            if ii==1 %  Forward finite-difference method
                Mdy_dz_vxyz(i,ii,k)= (-Mdz_vxyz(i,ii+2,k)+4.*...
                    Mdz_vxyz(i,ii+1,k)-3.*Mdz_vxyz(i,ii,k))/(2*hy); 
            elseif ii==2
                Mdy_dz_vxyz(i,ii,k)= (-Mdz_vxyz(i,ii+2,k)+4.*...
                    Mdz_vxyz(i,ii+1,k)-3.*Mdz_vxyz(i,ii,k))/(2*hy); 
            elseif ii==length(X)% Backward finite-difference method
                Mdy_dz_vxyz(i,ii,k)= (3.*Mdz_vxyz(i,ii,k)-4.*...
                    Mdz_vxyz(i,ii-1,k)+Mdz_vxyz(i,ii-2,k))/(2*hy);
            elseif ii==(length(X)-1)
                Mdy_dz_vxyz(i,ii,k)= (3.*Mdz_vxyz(i,ii,k)-4.*...
                    Mdz_vxyz(i,ii-1,k)+Mdz_vxyz(i,ii-2,k))/(2*hy);
            else % Centered finite-difference method
                Mdy_dz_vxyz(i,ii,k)= (-Mdz_vxyz(i,ii+2,k)+8.*...
                    Mdz_vxyz(i,ii+1,k)-8.*Mdz_vxyz(i,ii-1,k)+...
                    Mdz_vxyz(i,ii-2,k))/(12*hy);                                
            end
        end
    end
end

%--------- dz_dz_vxy calculation ---------
hz=diff(imersoes_2);
for i = 1:length(depth);
    for k = 1:length(X);
        for ii = 1:length(Y);
            if ii==1 %  Forward finite-difference method
                Mdz_dz_vxyz(i,k,ii)= (-Mdz_vxyz(i,k,ii+2)+4.*...
                    Mdz_vxyz(i,k,ii+1)-3.*Mdz_vxyz(i,k,ii))/...
                    (2*hz(i,1)); 
            elseif ii==2
                Mdz_dz_vxyz(i,k,ii)= (-Mdz_vxyz(i,k,ii+2)+4.*...
                    Mdz_vxyz(i,k,ii+1)-3.*Mdz_vxyz(i,k,ii))/...
                    (2*hz(i,1));
            elseif ii==length(Y)% Backward finite-difference method
                Mdz_dz_vxyz(i,k,ii)= (3.*Mdz_vxyz(i,k,ii)-4.*...
                    Mdz_vxyz(i,k,ii-1)+Mdz_vxyz(i,k,ii-2))/...
                    (2*hz(i,1));
            elseif ii==(length(Y)-1)
                Mdz_dz_vxyz(i,k,ii)= (3.*Mdz_vxyz(i,k,ii)-4.*...
                    Mdz_vxyz(i,k,ii-1)+Mdz_vxyz(i,k,ii-2))/...
                    (2*hz(i,1));
            else % Centered finite-difference method
                Mdz_dz_vxyz(i,k,ii)= (-Mdz_vxyz(i,k,ii+2)+8.*...
                    Mdz_vxyz(i,k,ii+1)-8.*Mdz_vxyz(i,k,ii-1)+...
                    Mdz_vxyz(i,k,ii-2))/(12*hz(i,1));                                
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Priting the graphics
disp('1 to print')
disp('0 to continue')
graphics=input('Do you wish to print the graphics: ');

if graphics ==1
    clear colorbar
    figure % % Bathymetry Map
    contourf(long_area,lat_area,prof_area(1:size(long_area(:)),...
        1:size(lat_area(:)))')
    hold on 
    scatter(posicao_geo3D(:,2),posicao_geo3D(:,1),'r','.','LineWidth',2);
    title ('Bathymetric map of the simulation area')
    xlabel ('Longitude [º]')
    ylabel ('Latitude [º]')
    colorbar=colorbar;
    colorbar.Label.String = 'Altitude in meters';
    set(gca,'Fontsize',10,'FontWeight','bold','LineWidth',2.0,'box','on');
    
    figure% 3D Graphic
    mesh(X_m(:,1),Y_m(1,:),prof_area(1:length(X_m),1:length(Y_m))') 
    title('3D Bathymetry')
    xlabel('Distance to the initial geographic point [m] - Longitude'); 
    ylabel('Distance to the initial geographic point [m] - Latitude'); 
    zlabel('Depth [m]');
    set(gca,'Fontsize',10,'FontWeight','bold','LineWidth',2.0,'box','on');
    
    figure
    scatter(posicao_geoS(:,2),posicao_geoS(:,1),'b');
    hold on 
    scatter(posicao_geoT(:,2),posicao_geoT(:,1),'b');
    hold on 
    scatter(posicao_geo3D(:,2),posicao_geo3D(:,1),'r','.','LineWidth',2);
    hold on
    scatter(posicao_geoS(Idx_S,2),posicao_geoS(Idx_S,1),'g',...
        's','LineWidth',2)
    hold on
    scatter(posicao_geoT(Idx_T,2),posicao_geoT(Idx_T,1),'+',...
        'LineWidth',2)
    title('Map with the distribution of points with the climatological data')
    xlabel ('Longitude [º]')
    ylabel ('Latitude [º]')
    set(gca,'Fontsize',10,'FontWeight','bold','LineWidth',2.0,'box','on');
elseif graphics ==0
    disp('End')
elseif isempty(graphics)==1
    disp('1 to print')
    disp('0 to continue')
    graphics=input('Do you wish to print the graphics: ');
else
    disp('1 to print')
    disp('0 to continue')
    graphics=input('Do you wish to print the graphics: ');
end