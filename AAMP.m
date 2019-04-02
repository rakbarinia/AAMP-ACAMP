function [mindist, minind] = MP2(X,m)%
    % X: time series, m: subsequence size
    % mindist: matrix profile, minind: the index of the nearest neighbor of each subsequence 
    [dim, Nb]=size(X); 
    s=Nb-m;
    Dmin=realmax*ones(1,s);
    minind=ones(1,s);

    for k=1:s-1
        %
        D=sum((X(1:m)-X(k+1:m+k)).^2);
        if D < Dmin(1)%1k
            Dmin(1)=D;
            minind(1)=k;
        end
        
        if D < Dmin(k)
            Dmin(k)=D;
            minind(k)=1;
        end

        for i=1:s-k%k=j-i
            kplusi=k+i;
            D=D-(X(i)-X(kplusi))^2 +(X(m+i)-X(m+kplusi))^2;
            
            if Dmin(i)> D
                minind(i)=kplusi;
                Dmin(i)=D;
            end
            
            if Dmin(kplusi)> D
                minind(kplusi)=i;
                Dmin(kplusi)=D;
            end
        end
  
    end
    mindist= sqrt(Dmin);
end
