function [mindist, minind] = MPnorme(X,m,p)%
    % X: time series, m: subsequence size
    % mindist: matrix profile, minind: the index of the nearest neighbor of each subsequence 
    [dim, Nb]=size(X); 
    s=Nb-m;
    Dmin=realmax*ones(1,s);
    minind=ones(1,s);

    for k=1:s-1
        D=sum((abs(X(1:m)-X(k+1:m+k))).^p);
        if D < Dmin(1)
            Dmin(1)=D;
            minind(1)=k;
        end
        
        if D < Dmin(k)
            Dmin(k)=D;
            minind(k)=1;
        end

        for i=1:s-k
            kplusi=k+i;
            D=D-(abs(X(i)-X(kplusi)))^p +(abs(X(m+i)-X(m+kplusi)))^p;
            
            if Dmin(i)> D
                minind(i)=kplusi;
                Dmin(i)=D;
            end
            %I use the symmetry of D to avoid doing a second loop.
            if Dmin(kplusi)> D
                minind(kplusi)=i;
                Dmin(kplusi)=D;
            end
        end
  
    end
    mindist=max(Dmin,0).^(1/p);
end
