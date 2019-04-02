function [mindist, minind] = MP2(X,m)%
    %X=serie temps, m la taille de la fenetre mobile, dist renvoie la sortie de
    %MP la plus petite distance (ie inf_j D_ij) et ind l indice associé
    [dim, Nb]=size(X); %dim la dimension, dim=1 classiquement et pour l instant, Nb= taille globale de la série
    s=Nb-m;%=taille de la fonction MP attention requiere d avoir Nb>m
    Dmin=realmax*ones(1,s);%mindist=realmax*ones(1,s);%initialisation de la distance (utile?Oui, inf à la place de realmax?)
    minind=ones(1,s);

    for k=1:s-1
        %on va calculer D_ij, on pose k=j-i, avt la deuxieme boucle i=1
        D=sum((X(1:m)-X(k+1:m+k)).^2);
        if D < Dmin(1)%1k
            Dmin(1)=D;
            minind(1)=k;
        end
        %symétrieFFmin(k)FFmin(k)FF
        if D < Dmin(k)
            Dmin(k)=D;
            minind(k)=1;
        end

        for i=1:s-k%k=j-i
            kplusi=k+i;
            D=D-(X(i)-X(kplusi))^2 +(X(m+i)-X(m+kplusi))^2;
            %Diiplusk=sqrt(2*m+2*F);
            if Dmin(i)> D
                minind(i)=kplusi;
                Dmin(i)=D;
            end
            %j utilise la symmétrie de D; Dij=Dji pour eviter de faire une
            %deuxième boucle
            if Dmin(kplusi)> D
                minind(kplusi)=i;
                Dmin(kplusi)=D;
            end
        end
  
    end
    mindist= sqrt(Dmin);%mindist= sqrt(2*m+2*s.*sqrt(abs(FFmin)));% ou ? a tester plus de fois;
end