function XGnm = XLaguerreGaussBeam(nu, mu, initialWaist, wavelength, units, r, th, z) 

    LPz = LaguerreParameters(z, initialWaist, wavelength, nu, mu, units);

    k    = LPz.k;
    Rz   = LPz.radius;
    wz   = LPz.waist;
    phiz = (2*nu+mu+1)*LPz.GouyPhase;
%     phiz = LPz.GouyPhase;

    argumentL = (2*r.^2)./(wz.^2);

    XGnm = XLaguerreG(nu,abs(mu),argumentL) ...
         .*exp( 1i*(k*r.^2./(2*Rz)))...
         .*exp(-1i*(phiz))...
         .*exp( 1i*mu*th)...
         .*exp( 1i*k*z)...
         ./(wz.^(mu+1));

    is_truncing = true;

    if is_truncing
        XGnm = XGnm.*exp(-(r/(initialWaist*sqrt(2*nu+mu+1))).^50);
    end

    XGnm(isnan(XGnm))=0;

end