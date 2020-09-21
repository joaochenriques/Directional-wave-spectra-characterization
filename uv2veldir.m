%converte velocidade leste(u) e velocidade norte (v) em magnitude (mag) e
%dire��o (dir) em coordenadas n�uticas
%[mag,dir]=uv2veldir(u,v) 
function [vel,dir]=uv2veldir(u,v)
vel=sqrt(u.^2+v.^2);
dir=rad2deg(atan2(u,v));    %dire��o n�o declinada em coordenadas n�uticas
dir(dir<0)=dir(dir<0)+360;  %transforma valores negativos em positivos
end
