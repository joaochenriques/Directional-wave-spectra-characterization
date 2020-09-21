%converte velocidade leste(u) e velocidade norte (v) em magnitude (mag) e
%direção (dir) em coordenadas náuticas
%[mag,dir]=uv2veldir(u,v) 
function [vel,dir]=uv2veldir(u,v)
vel=sqrt(u.^2+v.^2);
dir=rad2deg(atan2(u,v));    %direção não declinada em coordenadas náuticas
dir(dir<0)=dir(dir<0)+360;  %transforma valores negativos em positivos
end
