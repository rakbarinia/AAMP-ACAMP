function [mindist, minind] = MPz(X,m)%
    %X=serie temps, m la taille de la fenetre mobile, dist renvoie la sortie de
    %MP la plus petite distance (ie inf_j D_ij) et ind l indice associé
    [dim, Nb]=size(X); %dim la dimension, dim=1 classiquement et pour l instant, Nb= taille globale de la série
    s=Nb-m;%=taille de la fonction MP attention requiere d avoir Nb>m
    FFmin=realmax*ones(1,s);%mindist=realmax*ones(1,s);%initialisation de la distance (utile?Oui, inf à la place de realmax?)
    minind=ones(1,s);
    mm=1/m;
    A1=sum(X(1:m));%correspont  à la somme des elmts de la premiere fenetre (ie fenetre i); a changer pour la dim d idem pour les autres
    A21=sum(X(1:m).*X(1:m));%X(1:m)*X(1:m)';%sum(X(1:m).*X(1:m));%norm(X(1:m))^2; %sum(X(1:m).*X(1:m)); %avec des carres    
    x1=X(1);
    xm=X(1+m);
    B21=A21-x1^2+xm^2;    
    B1=A1-x1+xm;
    for k=1:s-1
        %on va calculer D_ij, on pose k=j-i, avt la deuxieme boucle i=1
        A=A1;%correspont  à la somme des elmts de la premiere fenetre (ie fenetre i); a changer pour la dim d idem pour les autres
        B=B1;%sum(X(k+1:k+m)); %correspont  à la somme des elmts de la deux fenetre (ie fenetre j)
        A2=A21;%X(1:m)*X(1:m)';%sum(X(1:m).*X(1:m));%norm(X(1:m))^2; %sum(X(1:m).*X(1:m)); %avec des carres
        B2=B21;%sum(X(k+1:k+m).*X(k+1:k+m));X(k+1:k+m)*X(k+1:k+m)';%norm(X(k+1:k+m));%sum(X(k+1:k+m).*X(k+1:k+m)); %idem
        kplus1=1+k;
        C=sum(X(kplus1:k+m).*X(1:m)); %somme des produit
        Z=A*B-m*C;
        FF=abs(Z)*(Z) /((A2-A^2*mm)*(B2-B^2*mm));
        %D1k=sqrt(2*m+2*F);;
        if FF < FFmin(1)%1k
            FFmin(1)=FF;
            minind(1)=k;
        end
        %symétrieFFmin(k)FFmin(k)FF
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
        %Diiplusk=sqrt(2*m+2*F);
        if FFmin(1)> FF
            minind(1)=kplus1;
            FFmin(1)=FF;
        end
        %j utilise la symmétrie de D; Dij=Dji pour eviter de faire une
        %deuxième boucle
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
            %Diiplusk=sqrt(2*m+2*F);
            if FFmin(i)> FF
                minind(i)=iplusk;
                FFmin(i)=FF;
            end
            %j utilise la symmétrie de D; Dij=Dji pour eviter de faire une
            %deuxième boucle
            if FFmin(iplusk)> FF
                minind(iplusk)=i;
                FFmin(iplusk)=FF;
            end
        end
  
    end
    s=sign(FFmin);
    mindist= sqrt(2*m+2*s.*sqrt(s.*FFmin));%mindist= sqrt(2*m+2*s.*sqrt(abs(FFmin)));% ou ? a tester plus de fois;
end