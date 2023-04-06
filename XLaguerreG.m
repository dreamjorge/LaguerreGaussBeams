function [Xn]= XLaguerreG(n,m,x)

Xn=(-1)^(n+1)/((n+(m+1)/2)^(m/2))*exp(-x/2).*x.^(m/2).*XLaguerreAssociated(n,m,x);

end