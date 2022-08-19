function [Rn,i,r] = sup_reflection3D(normal_R_sup,v,ctemp2)
%===================================================================%
% This function calculates the reflection of the acoustic ray 
% in the ocean surface.
%  Rn = normal vector to the plane at the point of reflection (Matrix)
%  i = incident acoustic ray 
%  r = reflected acoustic ray
%===================================================================%
Rn=[1-2*(normal_R_sup(1,1)^2),-2*normal_R_sup(1,1)*...
    normal_R_sup(2,1),-2*normal_R_sup(1,1)*normal_R_sup(3,1);...
    -2*normal_R_sup(1,1)*normal_R_sup(2,1),1-2*...
    (normal_R_sup(2,1)^2),-2*normal_R_sup(2,1)*normal_R_sup(3,1);...
    -2*normal_R_sup(1,1)*normal_R_sup(3,1),-2*normal_R_sup(2,1)*...
    normal_R_sup(3,1),1-2*(normal_R_sup(3,1)^2)];
Rn=Rn/norm(Rn);

i=[v(4,ctemp2);v(5,ctemp2);v(6,ctemp2)];
i=i/norm(i);

r=Rn*i;
r=r/norm(r);

end

