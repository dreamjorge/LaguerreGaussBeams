function [Hr,Hth,Hz] = HankelLaguerreGaussrthz(nu, mu, initialWaist, wavelength, ...
                                               r, th, z, ...
                                               ri, thi, zi, ...
                                               units, ... 
                                               hankeltype) 

    [Lr]  = LaguerreGaussBeam(nu, mu, initialWaist, wavelength, units, r,  thi, zi);
    [Lth] = LaguerreGaussBeam(nu, mu, initialWaist, wavelength, units, ri, th,  zi);
    [Lz]  = LaguerreGaussBeam(nu, mu, initialWaist, wavelength, units, ri, thi, z);

    [Xr]  = XLaguerreGaussBeam(nu, mu, initialWaist, wavelength, units, r,  thi, zi);
    [Xth] = XLaguerreGaussBeam(nu, mu, initialWaist, wavelength, units, ri, th,  zi);
    [Xz]  = XLaguerreGaussBeam(nu, mu, initialWaist, wavelength, units, ri, thi, z);

    if hankeltype == 1
        signhankel = +1;
    else
        signhankel = -1;
    end

    Hr  = Lr  + signhankel*1i*Xr;
    Hth = Lth + signhankel*1i*Xth;
    Hz  = Lz  + signhankel*1i*Xz; 

end