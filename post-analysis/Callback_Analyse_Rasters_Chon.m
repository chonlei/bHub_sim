
%Raster: input raw raster
test = 3;
if test==1
    load('temp_raster');
    Raster = a;
elseif test==2
    Raster = csvread('recording1_Raster.csv');
elseif test==3
    Raster = Ca_bi';
end

Signaux = struct;
Images = struct;
%Images.Nombre_Regions_Interet = number of cells (ROI)
Images.Nombre_Regions_Interet = size(Raster,2);
Signaux.Test_Total = ones(1,Images.Nombre_Regions_Interet);

Raster_Brasse=zeros(size(Raster));
size(Raster)
% Pour chaque zone d'interet fabrication d'un signal brasse
for i=1:Images.Nombre_Regions_Interet
    if(Signaux.Test_Total(1,i)==1)
        Activite_Tmp=Raster(:,i)';
        Bornes_Activite=Activite_Tmp(2:end)-Activite_Tmp(1:(end-1));
        if Activite_Tmp(1,1)==1
            Bornes_Activite=[ 1 Bornes_Activite];
        else
            Bornes_Activite=[0 Bornes_Activite];
        end
        if Activite_Tmp(1,end)==1 && Bornes_Activite(1,end)~=1
            Bornes_Activite(1,end)=-1;
        end
        if Activite_Tmp(1,end)==1 && Bornes_Activite(1,end)==1
            Bornes_Activite(1,end)=0;
        end

        Debut_Activite=find(Bornes_Activite==1);
        Fin_Activite=find(Bornes_Activite==-1);
        Longueur_Activite=Fin_Activite-Debut_Activite;
        Longueur_Inactivite=[Debut_Activite(1,1) Debut_Activite(2:end)-Fin_Activite(1:(end-1))];
        if Activite_Tmp(1,end)==0
            Longueur_Inactivite=[Debut_Activite(1,1) Debut_Activite(2:end)-Fin_Activite(1:(end-1)) size(Activite_Tmp,2)-Fin_Activite(1,end)];
        end
        
        Longueur_Activite_perm=Longueur_Activite(randperm(size(Longueur_Activite,2)));
        Longueur_Inactivite_perm=Longueur_Inactivite(randperm(size(Longueur_Inactivite,2)));
        index=1;
        for k=1:(size(Longueur_Inactivite,2)-1)
                Raster_Brasse(index:(index+Longueur_Inactivite_perm(1,k)-1),i)=0;
                Raster_Brasse((index+Longueur_Inactivite_perm(1,k)):(index+Longueur_Inactivite_perm(1,k)+Longueur_Activite_perm(1,k)-1),i)=1;
                index=index+Longueur_Inactivite_perm(1,k)+Longueur_Activite_perm(1,k);
        end
        Nombre_Evenements(1,i)=size(Debut_Activite,2);
% 
%         Debut_Activite_Brassee=floor(rand(1,Nombre_Evenements(1,i))*size(Raster,1))+1;
%         for j=1:Nombre_Evenements(1,i)
%             Raster_Brasse(Debut_Activite_Brassee(1,j):min(Debut_Activite_Brassee(1,j)+Longueur_Activite(1,j),size(Raster,1)),i)=1;
%         end
    end
end
    

% Boucle sur toutes les paires de cellules pour déterminer un indice de
% corrélation
Signaux.Correlation=zeros(Images.Nombre_Regions_Interet,Images.Nombre_Regions_Interet);
for i=1:Images.Nombre_Regions_Interet
    Signaux.Correlation(i,i)=1;
    for j=(i+1):Images.Nombre_Regions_Interet
        if ((Signaux.Test_Total(1,i)==1)&&(Signaux.Test_Total(1,j)==1))
            Synchro_Tmp=(Raster(:,i).*Raster(:,j))';
            Nombre_Synchro=size(find(Synchro_Tmp(2:end)-Synchro_Tmp(1:(end-1))==1),2);
            if(Synchro_Tmp(1,1)==1)
                Nombre_Synchro=Nombre_Synchro+1;
            end
            Signaux.Correlation(i,j)=Nombre_Synchro/sqrt(Nombre_Evenements(1,i)*Nombre_Evenements(1,j));
            Signaux.Correlation(j,i)=Signaux.Correlation(i,j);
        end
    end
end

% Boucle pour toutes les pares de cellules pour determiner la correlation
% entre signaux brasses
Signaux.Correlation_Brassee=zeros(Images.Nombre_Regions_Interet,Images.Nombre_Regions_Interet);
for i=1:Images.Nombre_Regions_Interet
    Signaux.Correlation(i,i)=1;
    for j=(i+1):Images.Nombre_Regions_Interet
        if ((Signaux.Test_Total(1,i)==1)&&(Signaux.Test_Total(1,j)==1))
            Synchro_Tmp=(Raster_Brasse(:,i).*Raster_Brasse(:,j))';
            Nombre_Synchro=size(find(Synchro_Tmp(2:end)-Synchro_Tmp(1:(end-1))==1),2);
            if(Synchro_Tmp(1,1)==1)
                Nombre_Synchro=Nombre_Synchro+1;
            end
            Signaux.Correlation_Brassee(i,j)=Nombre_Synchro/sqrt(Nombre_Evenements(1,i)*Nombre_Evenements(1,j));
            Signaux.Correlation_Brassee(j,i)=Signaux.Correlation_Brassee(i,j);
        end
    end
end

% Determination de la matrice de correlation
Correlation_Aleatoire_Mean_Tmp=mean(mean(Signaux.Correlation_Brassee));
Correlation_Aleatoire_Std_Tmp=mean(std(Signaux.Correlation_Brassee));
Seuil_Significativite_Correlation=Correlation_Aleatoire_Mean_Tmp+2*Correlation_Aleatoire_Std_Tmp;

Signaux.Correlation_Test=Signaux.Correlation>Seuil_Significativite_Correlation;

%% quick plot
a = sum(Signaux.Correlation_Test,1);
a = a/max(a)*100;
[N,edges] = histcounts(a,20);
N = N/sum(N)*100;
x=(edges(1:end-1)+edges(2:end))/2;
loglog(x,N,'o')
