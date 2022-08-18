%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------%
%                             REOART 2D                             %
%-------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REOART 2D (Real Environment Ocean-Acoustic Ray Tracing) is an 
% ocean-acoustic ray tracing model with the ability to study sound 
% propagation in a real ocean environment. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this script is where the sound speed and bathymetry are calculated 
%-------------------------------------------------------------------%
clearvars -global; close all; clc;
global Mvxy Mdx_vxy Mdy_vxy Mdx_dx_vxy Mdy_dx_vxy Mdy_dy_vxy...
    Mdx_dy_vxy dis_m depth

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------Initial geographic position----------------------
lat_pos=39;
long_pos=-10.1;
azimute = 270;
distance = 43.1965443;
inc_dis = 0.5; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------%
%                         GEBCO DATABASE                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Selecting the area of the world. In this case is just the 
% North Atlantic area (gebco=1).
gebco = 1;  
if gebco==1
    longitude=ncread('GEBCO_AtlanticoNorte.nc','lon');
    latitude=ncread('GEBCO_AtlanticoNorte.nc','lat');
    elevation=ncread('GEBCO_AtlanticoNorte.nc','elevation');
end

area=5; % Change the number to change the simulation area 
if area==1
    lat_S=34.5; lat_N=44.5; long_W=-15.5; long_E=-7.5; % North 
elseif area==2                                         % Atlantic
    lat_S=30; lat_N=45; long_W=-40; long_E=-5; % Northeast Atlantic
elseif area==3
    lat_S=28.5; lat_N=36.5; long_W=-20.5; long_E=-13.5; % Madeira
elseif area==4
    lat_S=36; lat_N=40; long_W=-33; long_E=-22; % Açores
elseif area==5
    lat_S=34; lat_N=44; long_W=-15; long_E=-6; % Continental Portugal
elseif area==6
     lat_S=30; lat_N=43; long_W=-32; long_E=-6; % Portugal ZEE
elseif area==7
    lat_S=46; lat_N=61; long_W=-25; long_E=0; % Scotland
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------%
%                           GEODESIC LINE                           %
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

% Posição geográfica com incrementos
posicao_geo = [y_lat;x_long];
posicao_geo=posicao_geo';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------%
%                          TEOS-10 FUNCTIONS                        %
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

% Finding the nearest neighbors between data and geographic positions 
Idx_S = knnsearch(posicao_geoS(:,1:2),posicao_geo(:,1:2));
Idx_T = knnsearch(posicao_geoT(:,1:2),posicao_geo(:,1:2));

% Construction of Salinity and Temperature matrices for the geographic 
% position
for i=1:length(posicao_geo)
    Practical_Salinity(i,:)=S(Idx_S(i),:);
    Insitu_Temperature(i,:)=T(Idx_T(i),:);
end

% Depth values in meters
depth=[0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,...
    95,100,125,150,175,200,225,250,275,300,325,350,375,400,425,...
    450,475,500,550,600,650,700,750,800,850,900,950,1000,1050,...
    1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1600,1650,...
    1700,1750,1800,1850,1900,1950,2000,2100,2200,2300,2400,2500,...
    2600,2700,2800,2900,3000,3100,3200,3300,3400,3500,3600,3700,...
    3800,3900,4000,4100,4200,4300,4400,4500,4600,4700,4800,4900,...
    5000,5100,5200,5300,5400,5500];

% Pressure calculation
for i=1:length(posicao_geo)
     p(i,1:length(depth)) = gsw_p_from_z(-depth,posicao_geo(i,1));
end

% Absolute salinity calculation
SA = gsw_SA_from_SP(Practical_Salinity,p,posicao_geo(:,2),...
    posicao_geo(:,1));

% Conservative temperature calculation
CT = gsw_CT_from_t(SA,Insitu_Temperature,p);

% Sound speed calculation
sound_speed = gsw_sound_speed(SA,CT,p);

% Sound speed matrix
Mvxy =sound_speed';

%-------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Calculating spatial gradients of the sound speed          %
%-------------------------------------------------------------------%
% The derivatives: dx_vxy dy_vxy dx_dx_vxy dy_dx_vxy dx_dy_vxy 
%                  dy_dy_vxy
% x - distance
% y - depth

%-----Centered finite-difference formulas----
% First Derivative: f'(xi)= -f(xi+2)+8f(xi+1)-8f(xi-1)+f(xi-2)/12h
% Second Derivative: f''(xi)= -f(xi+2)+16f(xi+1)-30f(xi)+16f(xi-1)-...
%                              f(xi-2)/12h^2
% h is the increment used between geographic positions

%-----Backward finite-difference formulas----
% First Derivative: f'(xi)= 3f(xi)-4f(xi-1)+f(xi-2)/2h
% Second Derivative: f''(xi)= 2f(xi)-5f(xi-1)+4f(xi-2)-f(xi-3)/h^2

%-----Forward finite-difference formulas-----
% First Derivative: f'(xi)= -f(xi+2)+4f(xi+1)-3f(xi)/2h
% Second Derivative: f''(xi)= -f(xi+3)+4f(xi+2)-5f(xi+1)+2f(xi)/h^2

%-------------------------------------------------------------------%
%%%%%%%%%%$%%%%%%%%%%%%% First Derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%
% dx - variation along the lines
% dy - variation along the columns

depth_2=[0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,...
    95,100,125,150,175,200,225,250,275,300,325,350,375,400,425,...
    450,475,500,550,600,650,700,750,800,850,900,950,1000,1050,...
    1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1600,1650,...
    1700,1750,1800,1850,1900,1950,2000,2100,2200,2300,2400,2500,...
    2600,2700,2800,2900,3000,3100,3200,3300,3400,3500,3600,3700,...
    3800,3900,4000,4100,4200,4300,4400,4500,4600,4700,4800,4900,...
    5000,5100,5200,5300,5400,5500,5600]; %acrecentei 5600 metros
depth_2=depth_2';

%--------- dx_vxy calculation (horizontal variation)----------
%h=diff(dis_m); 
hx=inc_dis*1852;
for i = 1:length(Mvxy);
    for ii = 1:length(posicao_geo);
        if ii==1 % Forward finite-difference method
            Mdx_vxy(i,ii)= (-Mvxy(i,ii+2)+4.*Mvxy(i,ii+1)-3.*...
                Mvxy(i,ii))/(2*hx); 
        elseif ii==2
            Mdx_vxy(i,ii)= (-Mvxy(i,ii+2)+4.*Mvxy(i,ii+1)-3.*...
                Mvxy(i,ii))/(2*hx); 
        elseif ii==length(posicao_geo)% Backward finite-difference method
            Mdx_vxy(i,ii)= (3.*Mvxy(i,ii)-4.*Mvxy(i,ii-1)+...
                Mvxy(i,ii-2))/(2*hx);
        elseif ii==(length(posicao_geo)-1)
            Mdx_vxy(i,ii)= (3.*Mvxy(i,ii)-4.*Mvxy(i,ii-1)+...
                Mvxy(i,ii-2))/(2*hx);
        else % Centered finite-difference method
            Mdx_vxy(i,ii)= (-Mvxy(i,ii+2)+8.*Mvxy(i,ii+1)-8.*...
                Mvxy(i,ii-1)+Mvxy(i,ii-2))/(12*hx);                                
        end
    end
   
end


%--------- dy_vxy calculation (vertical variation)------------
hy=diff(depth_2);
for i = 1:length(Mvxy);
    for ii = 1:length(posicao_geo);
        if i==1 % Forward finite-difference method
            Mdy_vxy(i,ii)= (-Mvxy(i+2,ii)+4.*Mvxy(i+1,ii)-3.*...
                Mvxy(i,ii))/(2*hy(i,1)); 
        elseif i==2
            Mdy_vxy(i,ii)= (-Mvxy(i+2,ii)+4.*Mvxy(i+1,ii)-3.*...
                Mvxy(i,ii))/(2*hy(i,1));
        elseif i==length(Mvxy)% Backward finite-difference method
            Mdy_vxy(i,ii)= (3.*Mvxy(i,ii)-4.*Mvxy(i-1,ii)+...
                Mvxy(i-2,ii))/(2*hy(i,1));
        elseif i==(length(Mvxy)-1)
            Mdy_vxy(i,ii)= (3.*Mvxy(i,ii)-4.*Mvxy(i-1,ii)+...
            Mvxy(i-2,ii))/(2*hy(i,1));
        else % Centered finite-difference method
            Mdy_vxy(i,ii)= (-Mvxy(i+2,ii)+8.*Mvxy(i+1,ii)-8.*...
                Mvxy(i-1,ii)+Mvxy(i-2,ii))/(12*hy(i,1));                                
        end
    end
   
end


%-------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%% Second Derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%
% The derivatives: dx_dx_vxy dy_dx_vxy dx_dy_vxy dy_dy_vxy
%-------------------------------------------------------------------%

%--------- dx_dx_vxy calculation --------- 
for i = 1:length(Mdx_vxy);
    for ii = 1:length(posicao_geo);
        if ii==1 % Forward finite-difference method
            Mdx_dx_vxy(i,ii)= (-Mdx_vxy(i,ii+2)+4.*Mdx_vxy(i,ii+1)...
                -3.*Mdx_vxy(i,ii))/(2*hx); 
        elseif ii==2
            Mdx_dx_vxy(i,ii)= (-Mdx_vxy(i,ii+2)+4.*Mdx_vxy(i,ii+1)...
                -3.*Mdx_vxy(i,ii))/(2*hx); 
        elseif ii==length(posicao_geo)% Backward finite-difference method
            Mdx_dx_vxy(i,ii)= (3.*Mdx_vxy(i,ii)-4.*Mdx_vxy(i,ii-1)+...
                Mdx_vxy(i,ii-2))/(2*hx);
        elseif ii==(length(posicao_geo)-1)
            Mdx_dx_vxy(i,ii)= (3.*Mdx_vxy(i,ii)-4.*Mdx_vxy(i,ii-1)+...
                Mdx_vxy(i,ii-2))/(2*hx);
        else % Centered finite-difference method
            Mdx_dx_vxy(i,ii)= (-Mdx_vxy(i,ii+2)+8.*Mdx_vxy(i,ii+1)...
                -8.*Mdx_vxy(i,ii-1)+Mdx_vxy(i,ii-2))/(12*hx);                                
        end
    end
   
end

%--------- dy_dx_vxy calculation ---------
for i = 1:length(Mdx_vxy);
    for ii = 1:length(posicao_geo);
        if i==1 % Forward finite-difference method
            Mdy_dx_vxy(i,ii)= (-Mdx_vxy(i+2,ii)+4.*Mdx_vxy(i+1,ii)...
                -3.*Mdx_vxy(i,ii))/(2*hy(i,1)); 
        elseif i==2
            Mdy_dx_vxy(i,ii)= (-Mdx_vxy(i+2,ii)+4.*Mdx_vxy(i+1,ii)...
                -3.*Mdx_vxy(i,ii))/(2*hy(i,1));
        elseif i==length(Mdx_vxy)% Backward finite-difference method
            Mdy_dx_vxy(i,ii)= (3.*Mdx_vxy(i,ii)-4.*Mdx_vxy(i-1,ii)+...
                Mdx_vxy(i-2,ii))/(2*hy(i,1));
        elseif i==(length(Mdx_vxy)-1)
            Mdy_dx_vxy(i,ii)= (3.*Mdx_vxy(i,ii)-4.*Mdx_vxy(i-1,ii)+...
                Mdx_vxy(i-2,ii))/(2*hy(i,1));
        else % Centered finite-difference method
            Mdy_dx_vxy(i,ii)= (-Mdx_vxy(i+2,ii)+8.*Mdx_vxy(i+1,ii)...
                -8.*Mdx_vxy(i-1,ii)+Mdx_vxy(i-2,ii))/(12*hy(i,1));                                
        end
    end
   
end

%--------- dx_dy_vxy calculation---------
for i = 1:length(Mdy_vxy);
    for ii = 1:length(posicao_geo);
        if ii==1 % Forward finite-difference method
            Mdx_dy_vxy(i,ii)= (-Mdy_vxy(i,ii+2)+4.*Mdy_vxy(i,ii+1)...
                -3.*Mdy_vxy(i,ii))/(2*hx); 
        elseif ii==2
            Mdx_dy_vxy(i,ii)= (-Mdy_vxy(i,ii+2)+4.*Mdy_vxy(i,ii+1)...
                -3.*Mdy_vxy(i,ii))/(2*hx); 
        elseif ii==length(posicao_geo)% Backward finite-difference method
            Mdx_dy_vxy(i,ii)= (3.*Mdy_vxy(i,ii)-4.*Mdy_vxy(i,ii-1)+...
                Mdy_vxy(i,ii-2))/(2*hx);
        elseif ii==(length(posicao_geo)-1)
            Mdx_dy_vxy(i,ii)= (3.*Mdy_vxy(i,ii)-4.*Mdy_vxy(i,ii-1)+...
                Mdy_vxy(i,ii-2))/(2*hx);
        else % Centered finite-difference method
            Mdx_dy_vxy(i,ii)= (-Mdy_vxy(i,ii+2)+8.*Mdy_vxy(i,ii+1)...
                -8.*Mdy_vxy(i,ii-1)+Mdy_vxy(i,ii-2))/(12*hx);                                
        end
    end
   
end

%--------- dy_dy_vxy calculation ---------
for i = 1:length(Mdy_vxy);
    for ii = 1:length(posicao_geo);
        if i==1 % Forward finite-difference method
            Mdy_dy_vxy(i,ii)= (-Mdy_vxy(i+2,ii)+4.*Mdy_vxy(i+1,ii)...
                -3.*Mdy_vxy(i,ii))/(2*hy(i,1)); 
        elseif i==2
            Mdy_dy_vxy(i,ii)= (-Mdy_vxy(i+2,ii)+4.*Mdy_vxy(i+1,ii)...
                -3.*Mdy_vxy(i,ii))/(2*hy(i,1));
        elseif i==length(Mdy_vxy)% Backward finite-difference method
            Mdy_dy_vxy(i,ii)= (3.*Mdy_vxy(i,ii)-4.*Mdy_vxy(i-1,ii)+...
                Mdy_vxy(i-2,ii))/(2*hy(i,1));
        elseif i==(length(Mdy_vxy)-1)
            Mdy_dy_vxy(i,ii)= (3.*Mdy_vxy(i,ii)-4.*Mdy_vxy(i-1,ii)+...
                Mdy_vxy(i-2,ii))/(2*hy(i,1));
        else % Centered finite-difference method
            Mdy_dy_vxy(i,ii)= (-Mdy_vxy(i+2,ii)+8.*Mdy_vxy(i+1,ii)...
                -8.*Mdy_vxy(i-1,ii)+Mdy_vxy(i-2,ii))/(12*hy(i,1));                                
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
    figure % Bathymetry Map
    contourf(long_area,lat_area,prof_area') 
    title ('Bathymetric map of the simulation area')
    xlabel ('Longitude [º]')
    ylabel ('Latitude [º]')
    hold on
    plot(posicao_geo(:,2),posicao_geo(:,1),'r.-')
    hold on 
    plot(posicao_geo(1,2),posicao_geo(1,1),'o-r')
    colorbar;
    colorbar=colorbar;
    colorbar.Label.String = 'Altitude in meters';
    set(gca,'Fontsize',10,'FontWeight','bold','LineWidth',2.0,'box','on')
    
    figure % 2D Bathymetry 
    plot(dis_km,perfil_fundo)
    hold on 
    plot(dis_km(1),perfil_fundo(1),'o-r')
    axis([0 max(dis_km) min(perfil_fundo) 0])
    xlabel('Range [km]')
    ylabel('Depth [m]')
    title(sprintf ...
     ('Bathymetry from lat. %0.1f and long. %0.1f on azimuth %0.0f', ...
     lat_pos,long_pos,azimute))
    set(gca,'Fontsize',10,'FontWeight','bold','LineWidth',2.0,'box','on');
    
    figure % 3D Bathymetry
    mesh(long_area,lat_area,prof_area') % Gráfico 3D
    title ('3D bathymetry map of the simulation area')
    xlabel ('Longitude [º]')
    ylabel ('Latitude [º]')
    set(gca,'Fontsize',10,'FontWeight','bold','LineWidth',2.0,'box','on')
    
    figure % Climatology - nearest neighbors
    scatter(posicao_geoS(:,2),posicao_geoS(:,1),'b');
    hold on 
    scatter(posicao_geoT(:,2),posicao_geoT(:,1),'b');
    hold on 
    scatter(posicao_geo(:,2),posicao_geo(:,1),'r','.','LineWidth',2);
    hold on
    scatter(posicao_geoS(Idx_S,2),posicao_geoS(Idx_S,1),'g','s',...
        'LineWidth',2)
    hold on
    scatter(posicao_geoT(Idx_T,2),posicao_geoT(Idx_T,1),'+','LineWidth',2)
    title('Map with the distribution of points with the climatological data')
    xlabel ('Longitude [º]')
    ylabel ('Latitude [º]') 
    set(gca,'Fontsize',10,'FontWeight','bold','LineWidth',2.0,'box','on')
    disp(' ')
    disp(' ')
    disp(' ')
    disp('Graphics printed')
    disp('End')
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