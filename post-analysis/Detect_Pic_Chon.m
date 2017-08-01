
Signaux = struct;

% Signal_Lisse_Input = Ca input
Signal_Lisse_Input = [0.99 1.2 4.1 12. 3 2 1 0 0.3];

Signal = Signal_Lisse_Input;

alpha=0.2;
beta=1;

%seuil_aires_pics=0.9;

Signaux.Activite_Courant=zeros(1,size(Signal,2));

Signal_Lisse_Tmp= Signal;
a=min(Signal_Lisse_Tmp);
b=mean(Signal_Lisse_Tmp);
c=max(Signal_Lisse_Tmp);

Signaux.Base_Signal_Lisse_Courant=mean(Signal_Lisse_Tmp(find(Signal_Lisse_Tmp<a+0.9*(b-a))))*ones(1,size(Signal_Lisse_Input,2));
%Signaux.Base_Signal_Lisse_Courant=mean(Signal_Lisse_Tmp(find(Signal_Lisse_Tmp<a+0.9*(b-a))))*ones(1,size(Signal_Lisse_Input,2));
%Signaux.Base_Signal_Lisse_Courant=Signal_Residu_Input-mean(Signal_Residu_Input)+Signaux.Base_Signal_Lisse_Courant;

%Test_Pic_Tmp=Signaux.Base_Signal_Lisse_Courant+alpha*std(Signal_Lisse_Tmp);
Test_Pic_Tmp=Signaux.Base_Signal_Lisse_Courant+alpha*(c-a);

Signaux.Activite_Totale_Courante=Signal_Lisse_Input>Test_Pic_Tmp;
    
