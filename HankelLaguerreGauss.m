function HLnm = HankelLaguerreGauss(nu, mu, initialWaist, wavelength, ...
                                    r, th, z, ...
                                    units, ... 
                                    hankeltype) 

    LGnm  =  LaguerreGaussBeam(nu, mu, initialWaist, wavelength, units, r,  th, z);
    XLGnm = XLaguerreGaussBeam(nu, mu, initialWaist, wavelength, units, r,  th, z);

    if hankeltype == 1
        signhankel = +1;
    else
        signhankel = -1;
    end

    HLnm = LGnm+1i*signhankel*XLGnm;

end