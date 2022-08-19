function [Rn,i,r] = sup_reflection2D(v,ctemp1)
%===================================================================%
% This function calculates the reflection of the acoustic ray in 
% the ocean surface
%  Rn = normal vector to the plane at the point of reflection(Matrix)
%  i = incident acoustic ray 
%  r = reflected acoustic ray
%===================================================================%

Rn=[-1,0;...
    0,-1];
Rn=Rn/norm(Rn);

i=[v(3,ctemp1);v(4,ctemp1)];
i=i/norm(i);

r=Rn*i;
r=r/norm(r);
end

