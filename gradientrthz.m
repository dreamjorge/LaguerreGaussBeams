function [mzr,mzth,mrth]=gradientrthz(fr,fth,fz,k,dr,dth,dz,r,th,z) 

    gz  = gradient(fz)/dz+k;
    gth = (1./r).*gradient(fth)./dth;
    gr  = gradient(fr)/dr;
    
    Nr  = size(gr,2);
    Nth = size(gth,2);
    Nz  = size(gz,2);
    
    index_r  = Nr/2 + 1 + floor(r/dr);
    index_th = Nth/2 + 1 + floor(th/dth);
    index_z  = floor(z/dz)+1;
    
    mzr  = gz(index_z)/gr(index_r);
    mzth = gz(index_z)/gth(index_th);
    mrth = gr(index_r)/gth(index_th);

end

