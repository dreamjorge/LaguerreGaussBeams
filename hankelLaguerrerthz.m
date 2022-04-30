function [H1r,H1th,H1z] = hankelLaguerrerthz(nu,mu,initialWaist,RayleighDistance,r,th,z,ri,thi,zi,hankeltype) 

    fthi  = exp(1i*mu*thi);

    fth   = exp(1i*mu*th);

    friz  = hankelLaguerrerz(nu,mu,initialWaist,RayleighDistance,ri,z , hankeltype);

    frizi = hankelLaguerrerz(nu,mu,initialWaist,RayleighDistance,ri,zi, hankeltype);

    frzi  = hankelLaguerrerz(nu,mu,initialWaist,RayleighDistance,r ,zi, hankeltype);
      
    H1r  = fthi.*frzi;

    H1th = fth.*frizi;

    H1z =  fthi.*friz;


end