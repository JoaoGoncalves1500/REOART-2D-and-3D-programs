function [normal_R] = normal_ponto_reflexao_2D(normals,n1)
%===================================================================%
% calculation of the normal vector at the surface at the 
% reflection point 
%===================================================================%

normal_R=[normals(n1(1,1),1)+normals(n1(1,2),1)/length(n1);...
    normals(n1(1,1),2)+normals(n1(1,2),2)/length(n1)];

normal_R=(normal_R/norm(normal_R));
normal_R(2,1)=-normal_R(2,1);
end

