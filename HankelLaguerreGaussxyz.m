function [Hx,Hy,Hz] = HankelLaguerreGaussxyz(nu, mu, initialWaist, wavelength, ...
                                             x, y, z, ...
                                             xi, yi, zi, ...
                                             units, ... 
                                             hankeltype) 
    [r_y, th_y] = cart2pol(xi,y);
    [r_x, th_x] = cart2pol(x,yi);
    [r_c, th_c] = cart2pol(xi,yi);


    Lx =  LaguerreGaussBeam(nu, mu, initialWaist, wavelength, units, r_x,  th_x, zi);
    Ly =  LaguerreGaussBeam(nu, mu, initialWaist, wavelength, units, r_y,  th_y, zi);
    Lz =  LaguerreGaussBeam(nu, mu, initialWaist, wavelength, units, r_c,  th_c, z);

    XLx = XLaguerreGaussBeam(nu, mu, initialWaist, wavelength, units, r_x,  th_x, zi);
    XLy = XLaguerreGaussBeam(nu, mu, initialWaist, wavelength, units, r_y,  th_y, zi);
    XLz = XLaguerreGaussBeam(nu, mu, initialWaist, wavelength, units, r_c,  th_c, z);


    if hankeltype == 1
        signhankel = +1;
    else
        signhankel = -1;
    end

    Hx  = Lx + signhankel*1i*XLx;
    Hy  = Ly + signhankel*1i*XLy;
    Hz  = Lz + signhankel*1i*XLz; 

end