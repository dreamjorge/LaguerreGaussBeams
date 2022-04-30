function frz = hankelLaguerrerz(nu,mu,initialWaist,RayleighDistance,r,zi, hankelType)


if hankelType==1

    frz=       laguerregz(nu,mu,initialWaist,RayleighDistance,r,zi)...
        + 1i*pxlaguerregz(nu,mu,initialWaist,RayleighDistance,r,zi);

elseif hankelType==2

    frz =       laguerregz(nu,mu,initialWaist,RayleighDistance,r,zi)...
        - 1i*pxlaguerregz(nu,mu,initialWaist,RayleighDistance,r,zi);

end