%converte velocidades magnitude e dire��o em u,v
%direcao em coordenadas n�uticas
%[u,v]=magdir2uv(mag,dir,declinacao)
function [u,v]=magdir2uv(mag,dir,declinacao)
dir=dir+declinacao;
u=mag.*sin(deg2rad(dir));
v=mag.*cos(deg2rad(dir));
end
