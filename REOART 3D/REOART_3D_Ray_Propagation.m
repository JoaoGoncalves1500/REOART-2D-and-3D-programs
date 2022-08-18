%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------%
%                    REOART 3D - Ray Propagation                    %
%-------------------------------------------------------------------%
% This is the REOART 3D main script where the propagation of the 
% acoutic ray is calculated
%-------------------------------------------------------------------%
%
clc
Vmatrix=[]; % ray trajectory matrix
Tsmatrix=[]; % propagation time matrix

% Initial conditions (change as you intend)
delta=0.5;
ts=0:delta:200;
theta0=0; %Angle in degrees 0<theta0<90
theta0=pi*theta0/180;
phi0=3; %Angle in degrees -90<phi0<90
phi0=pi*phi0/180;
%
a0= 50000;
b0= 80000;
c0=500;
% 
x0=cos(phi0)*cos(theta0);
y0=cos(phi0)*sin(theta0);
z0=sin(phi0);
%
v0 =[a0;b0;c0;x0;y0;z0];

Zz=-prof_area(1:length(X_m),1:length(Y_m))';

for il=1:20
    v=dae4('pa3dr',ts,v0,10);
    Zq=interp2(Y_m,X_m,Zz,v(1,:),v(2,:));
    ctemp1=min(find(v(3,2:end)>= Zq(1,2:end)));
    ctemp2=min(find(v(3,2:end)<=0));
    if isempty(ctemp1)==1 && isempty(ctemp2)==1
        disp('No Reflection')
        ts=ts+ts(end);
        Vmatrix=[Vmatrix v(:,1:end-1)];
        Tsmatrix=[Tsmatrix ts(1,1:end-1)];
        v0=v(:,end);
        
    elseif isempty(ctemp1)~=1 && isempty(ctemp2)==1
        disp('Bottom Reflection')
        %
        Zq_lat=interp2(Y_m,X_m,Zz,v(1,ctemp1+1)+100,v(2,ctemp1));
        Zq_long=interp2(Y_m,X_m,Zz,v(1,ctemp1),v(2,ctemp1+1)+100);
        %
        p1=[v(1,ctemp1),v(2,ctemp1),v(3,ctemp1)];
        p2=[v(1,ctemp1+1)+100,v(2,ctemp1),Zq_lat];
        p3=[v(1,ctemp1),v(2,ctemp1+1)+100,Zq_long];
        %
        A=p1-p2;
        B=p3-p2;
        C = cross(A,B);
        C = C./norm(C);
        dot(C,A);
        dot(C,B);
        normal_R=C';
        %
        [Rn,i,r]=reflection3D(normal_R,v,ctemp1) 
        %
        a0=v(1,ctemp1);
        b0=v(2,ctemp1);
        c0=v(3,ctemp1); 
        x0=r(1,1);
        y0=r(2,1);
        z0=r(3,1);
        %
        v0 =[a0;b0;c0;x0;y0;z0]; % New initial conditions
        Vmatrix=[Vmatrix v(:,1:ctemp1)];
        Tsmatrix=[Tsmatrix ts(1,1:ctemp1)];
        
    elseif isempty(ctemp1)==1 && isempty(ctemp2)~=1
        disp('Surface reflection')
        %
        p1=[v(1,ctemp2),v(2,ctemp2),v(3,ctemp2)];
        p2=[v(1,ctemp2+1)+100,v(2,ctemp2),v(3,ctemp2)];
        p3=[v(1,ctemp2),v(2,ctemp2+1)+100,v(3,ctemp2)];
        %
        A=p1-p2;
        B=p3-p2;
        C = cross(A,B);
        C = C./norm(C);
        dot(C,A)
        dot(C,B)
        normal_R_sup=C'
        %
        [Rn,i,r]=sup_reflection3D(normal_R_sup,v,ctemp2)
        %
        a0=v(1,ctemp2);
        b0=v(2,ctemp2);
        c0=v(3,ctemp2); 
        x0=r(1,1);
        y0=r(2,1);
        z0=r(3,1);
        %
        v0 =[a0;b0;c0;x0;y0;z0]; % New initial conditions
        Vmatrix=[Vmatrix v(:,1:ctemp2-1)];
        Tsmatrix=[Tsmatrix ts(1,1:ctemp2-1)];
    end 
    v0
end

figure
mesh(X_m(:,1),Y_m(1,:),prof_area(1:length(X_m),1:length(Y_m))')
% xlim([{number} {number}])        -| To see closely uncomment 
% ylim([{number} {number}])        -| the xlim and ylim. 
hold on;
plot3(Vmatrix(1,:),Vmatrix(2,:),-Vmatrix(3,:),'r','LineWidth',1.5)
set(gca,'Fontsize',10,'FontWeight','bold','LineWidth',2.0,'box','on');
a1=xlabel('Distance to the initial geographic point [m] - Longitude'); 
a2=ylabel('Distance to the initial geographic point [m] - Latitude'); 
a3=zlabel('Depth [m]');
a4=title('Acoustic ray propagation');
set(a1,'Fontsize',12,'FontWeight','bold');
set(a2,'Fontsize',12,'FontWeight','bold');
set(a3,'Fontsize',12,'FontWeight','bold');
set(a4,'Fontsize',12,'FontWeight','bold');
