% use solid (one) averaged threshold and baseline rather than individual threshold
function Ca_bi = binarise_trace2(Ca_input,savefile)
% Signal_Lisse_Input = Ca input
%Ca_input = dlmread('Ca_model_2_morphology_3_seed_1_mode_0_96x5001.txt');
Ca_bi = zeros(size(Ca_input));
Ca_max = mean(max(Ca_input));
Ca_min = mean(min(Ca_input));
Ca_mean = mean(mean(Ca_input));
MaxMinDiff = Ca_max - Ca_min;
baselineAll = mean(Ca_input(Ca_input<Ca_min+0.9*(Ca_mean-Ca_min)));
for cellid = 1:size(Ca_input,1)

Signaux = struct;

Signal_Lisse_Input = Ca_input(cellid,:);

Signal = Signal_Lisse_Input;

alpha=0.2;
beta=1;

%seuil_aires_pics=0.9;

Signaux.Activite_Courant=zeros(1,size(Signal,2));

Signal_Lisse_Tmp= Signal;
a=min(Signal_Lisse_Tmp);
b=mean(Signal_Lisse_Tmp);
c=max(Signal_Lisse_Tmp);

Signaux.Base_Signal_Lisse_Courant=baselineAll*ones(1,size(Signal_Lisse_Input,2));
%Signaux.Base_Signal_Lisse_Courant=mean(Signal_Lisse_Tmp(find(Signal_Lisse_Tmp<a+0.9*(b-a))))*ones(1,size(Signal_Lisse_Input,2));
%Signaux.Base_Signal_Lisse_Courant=Signal_Residu_Input-mean(Signal_Residu_Input)+Signaux.Base_Signal_Lisse_Courant;

%Test_Pic_Tmp=Signaux.Base_Signal_Lisse_Courant+alpha*std(Signal_Lisse_Tmp);
%Test_Pic_Tmp=Signaux.Base_Signal_Lisse_Courant+alpha*(c-a);
Test_Pic_Tmp=Signaux.Base_Signal_Lisse_Courant+alpha*MaxMinDiff;

Signaux.Activite_Totale_Courante=Signal_Lisse_Input>Test_Pic_Tmp;

Ca_bi(cellid,:) = Signaux.Activite_Totale_Courante;

end

if ~strcmp(savefile,'None')
    dlmwrite(savefile,Ca_bi,'\t')
end