%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------%
%                    REOART 2D - Ray Propagation                    %
%-------------------------------------------------------------------%
% This is the REOART 2D main script where the propagation of the 
% acoutic ray is calculated
%-------------------------------------------------------------------%
%
clc
Vmatrix=[]; % ray trajectory matrix
Tsmatrix=[]; % propagation time matrix

% Initial conditions (change as you intend)
delta=10;
ts=0:delta:1000*delta;
theta0=0; %Angle in degrees 0<theta0<90
theta0=pi*theta0/180;
x0=3500; % Horizontal position
y0=2500; % Vertical position 
z0=cos(theta0);
w0=sin(theta0);
v0=[x0;y0;z0;w0] % Ray initial conditions

for il=1:5
    v=dae4('pa2dr',ts,v0,10);
    Yq=interp1(dis_m,-perfil_fundo,v(1,:));
    ctemp1=min(find(v(2,2:end)<=1)); 
    ctemp2=min(find(v(2,2:end)>Yq(1,2:end)));
    
    if isempty(ctemp1)==1 && isempty(ctemp2)==1 % Without reflection
        disp('No reflection')
        ts=ts+ts(end);
        Vmatrix=[Vmatrix v(:,1:end-1)];
        Tsmatrix=[Tsmatrix ts(1,1:end-1)];
        v0=v(:,end);
        
    elseif isempty(ctemp1)~=1 && isempty(ctemp2)==1 % Surface reflection
        disp('Surface reflection')
        Vmatrix=[Vmatrix v(:,1:ctemp1)];
        Tsmatrix=[Tsmatrix ts(1,1:ctemp1)];
        v0=Vmatrix(:,end); % New initial conditions
        [Rn,i,r] = sup_reflection2D(v,ctemp1)
        v0(4)=r(2,1);
        ts=ts+Tsmatrix(end);
    elseif isempty(ctemp1)==1 && isempty(ctemp2)~=1 % Bottom reflection
        disp('Bottom reflection')
        Vmatrix=[Vmatrix v(:,1:ctemp2)];
        Tsmatrix=[Tsmatrix ts(1,1:ctemp2)];
        N=points'; 
        newpoint=[v(1,ctemp2),v(2,ctemp2)];
        [n,d] = knnsearch(N,newpoint,'k',3);  
        n1=unique(n);
        N_closest = N(n1,:); 
        [normal_R] = normal_ponto_reflexao_2D(normals',n1)
        [Rn,i,r] = reflection2D(normal_R,v,ctemp2)
        v0=Vmatrix(1:2,end); % New initial conditions
        v0(3,1)=r(1,1);
        v0(4,1)=r(2,1);
        ts=ts+Tsmatrix(end);  
    end
end

figure
plot(Vmatrix(1,:),Vmatrix(2,:),'r','LineWidth',1)
hold on
plot(dis_m,-perfil_fundo)
axis ij
grid on
set(gca,'Fontsize',10,'FontWeight','bold','LineWidth',2.0,'box','on');
a1=xlabel('Range [m]'); 
a2=ylabel('Depth [m]'); 
a3=title('Acoustic ray propagation');
set(a1,'Fontsize',12,'FontWeight','bold');
set(a2,'Fontsize',12,'FontWeight','bold')
set(a3,'Fontsize',12,'FontWeight','bold')

