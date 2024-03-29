% function [Ln]= LagG(n,m,x)
% 
% Ln=1./(gamma(m+1)*exp(x./2).*((n+(m+1)/2)*x).^(-m/2)).*hypergeom(-n,m+1,x);
% 
% end

function [Ln]= LaguerreG(n,m,x)

if n==0
    
    Ln=1;
    
else
    
    if n==1
        
        Ln=1+m-x;
        
    else
        
        l0=1;
        l1=1+m-x;
        
        for k=2:n
            
            Ln=((2*(k-1)+m+1-x).*l1-(k+m-1).*l0)./k;
            l0=l1;
            l1=Ln;
            
        end
        
    end
    
end

Cnm = (gamma(n+1)/gamma(n+m+1));
Ln = Cnm.*(exp(-x./2).*((n+(m+1)/2)*x).^(m/2)).*LaguerreAssociated(n,m,x);
            
end 
    
