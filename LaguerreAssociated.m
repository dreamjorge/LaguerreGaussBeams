function [Ln]= LaguerreAssociated(n,m,x)

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
