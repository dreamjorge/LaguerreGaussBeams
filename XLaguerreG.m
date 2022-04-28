function [Xn]= XLaguerreG(r,n,m,x)

Xn=exp(-x/2).*x.^(m/2).*XLaguerre(r,n,m,x);

end