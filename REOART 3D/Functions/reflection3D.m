function [Rn,i,r] = reflection3D(normal_R,v,ctemp1)
%===================================================================%
% This function calculates the reflection of the acoustic ray 
% in the bottom of the ocean.
%  Rn = normal vector to the plane at the point of reflection (Matrix)
%  i = incident acoustic ray 
%  r = reflected acoustic ray
%===================================================================%
Rn=[1-2*(normal_R(1,1)^2),-2*normal_R(1,1)*normal_R(2,1),-2*...
    normal_R(1,1)*normal_R(3,1);...
    -2*normal_R(1,1)*normal_R(2,1),1-2*(normal_R(2,1)^2),-2*...
    normal_R(2,1)*normal_R(3,1);...
    -2*normal_R(1,1)*normal_R(3,1),-2*normal_R(2,1)*...
    normal_R(3,1),1-2*(normal_R(3,1)^2)];
Rn=Rn/norm(Rn);

i=[v(4,ctemp1);v(5,ctemp1);v(6,ctemp1)];
i=i/norm(i);

r=Rn*i;
r=r/norm(r);
end

