%converte angulos de coordenadas cartesianas para coordenadas náuticas 
%out=cart2naut(in)
function out=cart2naut(in)
out=90-in;
out(out<0)=out(out<0)+360;
out(out>=360)=out(out>=360)-360;
