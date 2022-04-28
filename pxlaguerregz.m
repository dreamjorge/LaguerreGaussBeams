function [Xn]=xlaguerregz(nu,mu,wo,zl,x,z) 
    
k=2*zl/wo^2;
    Rz=z+zl^2./z;
    wz=wo*sqrt(1+(z./zl).^2);
    phiz=(2*nu+mu+1)*atan(z./zl);
%     Truncamiento ligado al plano z=0;

Xn=exp(-(x/(wo*sqrt(2*nu+mu+1))).^50).*XLaguerreG(97,nu,abs(mu),2*x.^2./(wz.^2)).*exp(1i*(k*x.^2./(2*Rz)-phiz));%./(wz.^(mu+1));
Xn(isnan(Xn))=0;
end