%converte velocidades magnitude e direção em u,v
%direcao em coordenadas náuticas
%[u,v]=magdir2uv(mag,dir,declinacao)
function [u,v]=magdir2uv(mag,dir,declinacao)
dir=dir+declinacao;
u=mag.*sin(deg2rad(dir));
v=mag.*cos(deg2rad(dir));
end
