function [mindist, minind] = MPz(X,m)%
    % X: time series, m: subsequence size
    % mindist: matrix profile, minind: the index of the nearest neighbor of each subsequence 
    [dim, Nb]=size(X); 
    s=Nb-m;
    FFmin=realmax*ones(1,s);
    minind=ones(1,s);
    mm=1/m;
    A1=sum(X(1:m));
    A21=sum(X(1:m).*X(1:m));    
    x1=X(1);
    xm=X(1+m);
    B21=A21-x1^2+xm^2;    
    B1=A1-x1+xm;
    for k=1:s-1
        A=A1;
        B=B1;
        A2=A21;
        B2=B21;
        kplus1=1+k;
        C=sum(X(kplus1:k+m).*X(1:m)); 
        Z=A*B-m*C;
        FF=abs(Z)*(Z) /((A2-A^2*mm)*(B2-B^2*mm));

        if FF < FFmin(1)%1k
            FFmin(1)=FF;
            minind(1)=k;
        end

        if FF < FFmin(k)
            FFmin(k)=FF;
            minind(k)=1;
        end
        A=A-x1 + xm;
        A2=A2-x1^2+xm^2;
        xk=X(kplus1);
        xkm=X(m+kplus1);
        B=B -xk + xkm;
        B2=B2 -xk^2 + xkm^2;
        B21=B2;
        B1=B;
        C= C-x1*xk+ xm*xkm;
        Z=A*B-m*C;
        FF=abs(Z)*(Z)/((A2-A^2*mm)*(B2-B^2*mm));
        
        if FFmin(1)> FF
            minind(1)=kplus1;
            FFmin(1)=FF;
        end
        
        if FFmin(kplus1)> FF
            minind(kplus1)=1;
            FFmin(kplus1)=FF;
        end
        for i=2:s-k%k=j-i
            xi=X(i);
            xmi=X(m+i);
            A=A-xi + xmi;
            A2=A2-xi^2 + xmi^2;
            iplusk=i+k;
            xik=X(iplusk);
            xmik=X(m+iplusk);
            B=B -xik+ xmik;
            B2=B2 -xik^2 + xmik^2;
            C= C-xi*xik+ xmi*xmik;
            Z=A*B-m*C;
            FF=abs(Z)*(Z)/((A2-A^2*mm)*(B2-B^2*mm));
            
            if FFmin(i)> FF
                minind(i)=iplusk;
                FFmin(i)=FF;
            end
            %I use the symmetry of D to avoid doing a second loop. 
            if FFmin(iplusk)> FF
                minind(iplusk)=i;
                FFmin(iplusk)=FF;
            end
        end
  
    end
    s=sign(FFmin);
    mindist= sqrt(2*m+2*s.*sqrt(s.*FFmin));
end
