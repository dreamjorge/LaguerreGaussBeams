function [mzr,mzth,mrth]=gradientrthz(fr,fth,fz,k,dr,dth,dz,r,th,z) 
%funcion que calcula el gradiente en un punto dado 
% f(x,y,z) a y,z constante en la primer componente
% f(x,y,z) a x,z constante en la segunda componente
% f(x,y,z) a x,y constante en la tercer componente

%derivadas parciales
gz=gradient(fz)/dz+k;
gth = (1./r).*gradient(fth)./dth;
gr=gradient(fr)/dr;
N=size(gr,2);
%pendientes
mzr  = gz(floor(z/dz)+1)/gr(N/2+1+floor(r/dr));
mzth = gz(floor(z/dz+1))    /gth(N/2+1+floor(th/dth));
mrth = gr(N/2+1+floor(r/dr))/gth(N/2+1+floor(th/dth));

