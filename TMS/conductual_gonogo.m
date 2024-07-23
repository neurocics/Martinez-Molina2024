clear 
clc

PP = '/Volumes/MariaPaz/DATOSENPROCESO/FONDECYT_R_2018';

Sujetos_A = {
'C_SL_19031984',...
'CADC_01031994',...
'CAMG_18121984',...
'CCAC_12071987',...
'CEMM_09041993',...
'CIPA_20031995',...
'CMST_12121996',...
'GBVO_25121991',...
'KKPM_16121998',...
'NELR_27121994',...
'NIFL_23091991',...
'TICA_05021995',...
'YTFV_09081990',...
'C_PB_09031993',...
'GADT_29111991',...
'NASO_18101996',...
'MFSG_02071992',...
'D_RM_31121985',...
'FECM_30051994',...
'GBMC_29041996',...
'MSGC_13081985',...
'CARO_15091995',...
'FJPL_05121998',...
}...
;
% 
% Sujetos_A_orden = {
% '1',...
% '1',...
% '1',...
% '1',...
% '1',...
% '1',...
% '1',...
% '1',...
% '1',...
% '1',...
% '1',...
% '1',...
% '1',...
% '1',...
% '2',...
% '2',...
% }...


load([  PP '/' Sujetos_A{1} '/LOG_A/RT_gonogo.mat' ] , 'RT')

RT.OTHER.TMS_SITE = cell(size(RT.OTHER.TMS));
RT.OTHER.TMS_SITE(:) = {'A'};
RT.OTHER.Sub(:) = Sujetos_B{1};

TRT = RT;


HITng50 = mean ( RT.rt( RT.est==21 & RT.OTHER.TMS==50 ) == -99  );
FAng50 = mean (RT.rt(RT.est==10 & RT.OTHER.TMS==50 ) == -99);
Dng50 = (exp(HITng50) - exp(FAng50));

HITng51 = mean ( RT.rt( RT.est==21 & RT.OTHER.TMS==51 ) == -99  );
FAng51 = mean (RT.rt(RT.est==10 & RT.OTHER.TMS==51 ) == -99);
Dng51 = (exp(HITng51) - exp(FAng51));

HITng52 = mean ( RT.rt( RT.est==21 & RT.OTHER.TMS==52 ) == -99  );
FAng52 = mean (RT.rt(RT.est==10 & RT.OTHER.TMS==52) == -99);
Dng52 = (exp(HITng52) - exp(FAng52));


HITng = mean (RT.rt(RT.est==21) == -99);
FAng = mean (RT.rt(RT.est==10) == -99);
Dng = (exp(HITng) - exp(FAng));

HITg = mean (RT.rt(RT.est==10) > 0);
FAg = mean (RT.rt(RT.est==21) > 0);
Dg = (exp(HITg) - exp(FAg));


clear RT






for s = 2:length(Sujetos_A)
    
    
load([  PP '/' Sujetos_A{s} '/LOG_A/RT_gonogo.mat' ] , 'RT')    

RT.OTHER.TMS_SITE = cell(size(RT.OTHER.TMS));
RT.OTHER.TMS_SITE(:) = {'A'}; 
RT.OTHER.Sub(:) = Sujetos_A{s};
TRT = rt_merge(RT, TRT,0); % El parametro 0 es para que no los desordene


HITng52 = mean ( RT.rt( RT.est==21 & RT.OTHER.TMS==52 ) == -99  );
FAng52 = mean (RT.rt(RT.est==10 & RT.OTHER.TMS==52) == -99);
Dng52(s) = (exp(HITng52) - exp(FAng52));


HITng51 = mean ( RT.rt( RT.est==21 & RT.OTHER.TMS==51 ) == -99  );
FAng51 = mean (RT.rt(RT.est==10 & RT.OTHER.TMS==51 ) == -99);
Dng51(s) = (exp(HITng51) - exp(FAng51));


clear RT
    
    
end

PP = '/Volumes/MariaPaz/DATOSENPROCESO/FONDECYT_R_2018';

Sujetos_B = {...
'C_SL_19031984',...
'CADC_01031994',...
'CAMG_18121984',...
'CCAC_12071987',...
'CEMM_09041993',...
'CIPA_20031995',...
....%'CMST_12121996',...
...%'GBVO_25121991',...
'KKPM_16121998',...
'NELR_27121994',...
'NIFL_23091991',...
'TICA_05021995',...
'YTFV_09081990',...
'C_PB_09031993',...
'GADT_29111991',...
'NASO_18101996',...
'MFSG_02071992',...
'D_RM_31121985',...
'FECM_30051994',...
'GBMC_29041996',...
'HJPF_13021991',...
'MSGC_13081985',...
'CARO_15091995',...
'FJPL_05121998',...
}...
;


% Sujetos_B_orden = {...
% '2',...
% '2',...
% '2',...
% '2',...
% '2',...
% '2',...
% ....%'CMST_12121996',...
% ...%'GBVO_25121991',...
% '2',...
% '2',...
% '2',...
% '2',...
% '2',...
% '2',...
% '1',...
% '1',...
% }...
% ;

for s = 1:length(Sujetos_B)
    
    
load([  PP '/' Sujetos_B{s} '/LOG_B/RT_gonogo.mat' ] , 'RT')    
    
RT.OTHER.TMS_SITE = cell(size(RT.OTHER.TMS));
RT.OTHER.TMS_SITE(:) = {'B'};
RT.OTHER.Sub(:) = Sujetos_B{s};


TRT = rt_merge(RT, TRT,0); % El parametro 0 es para que no los desordene


HITng52 = mean ( RT.rt( RT.est==21 & RT.OTHER.TMS==52 ) == -99  );
FAng52 = mean (RT.rt(RT.est==10 & RT.OTHER.TMS==52) == -99);
Dng52(s) = (exp(HITng52) - exp(FAng52));


HITng51 = mean ( RT.rt( RT.est==21 & RT.OTHER.TMS==51 ) == -99  );
FAng51 = mean (RT.rt(RT.est==10 & RT.OTHER.TMS==51 ) == -99);
Dng51(s) = (exp(HITng51) - exp(FAng51));


clear RT
    
    
end
%%

load('conductual_12122019.mat')


%%
TRT.OTHER.nnb=ones(size(TRT.OTHER.seq));
TRT.OTHER.nnb(:)=repmat(1:20,[1,numel(TRT.OTHER.seq)/20]);
%
Error  = (TRT.rt>10 & TRT.est==21) | (TRT.rt<0 & TRT.est==10);
pError = [ 0  Error(1:end-1) ] ;
Error = Error(TRT.rt>10 & TRT.est==10 )';
pError = pError(TRT.rt>10 & TRT.est==10 )';
Drt = TRT.rt( TRT.rt>10 & TRT.est==10 )';
Dseq = TRT.OTHER.seq( TRT.rt>10 & TRT.est==10 )';
Dsu = TRT.OTHER.Sub( TRT.rt>10 & TRT.est==10 )';
Dtms = TRT.OTHER.TMS( TRT.rt>10 & TRT.est==10 )';
DtmsB = ifcellis(  TRT.OTHER.TMS_SITE, 'B');
DtmsB = DtmsB( TRT.rt>10 & TRT.est==10 )';
DtmsA = ifcellis(  TRT.OTHER.TMS_SITE, 'A');
DtmsA = DtmsA( TRT.rt>10 & TRT.est==10 )';
Dtms51 = Dtms==51;
Dtms52 = Dtms==52;
Dtms50 = Dtms==50;
Dtms=Dtms51+Dtms52;


DATA = table(Drt,Dseq,Dsu,Dtms51,Dtms52,Dtms50,Dtms, DtmsB, DtmsA, Error, pError);



%%

Model = fitlme(DATA, 'Drt ~ Dseq + Dtms51 + Dtms50 + (1   |Dsu)  ')
% 98597    98640    -49293           98585   


Model = fitlme(DATA, 'Drt ~ Dseq  + Dtms51 + Dtms52 + (1 + Dtms51 + Dtms52 +   |Dsu)  ')
% 98597    98640    -49293           98585   



%%

Model = fitlme(DATA, 'Drt ~ Dseq + Dseq:Dtms + Dseq:Dtms51  + Dseq:Dtms:DtmsA + Dseq:Dtms51:DtmsA +  Dtms + Dtms51 + DtmsA + Dtms:DtmsA  + Dtms:Dtms51:DtmsA + (1 |Dsu)  ')
%98599    98642    -49294           98587   





% 
Model = fitlme(DATA, 'Drt ~ Dseq + Dseq:Dtms + Dseq:Dtms51  + Dseq:Dtms:DtmsA + Dseq:Dtms51:DtmsA+  Dtms + Dtms51 + DtmsA + Dtms:DtmsA  + Dtms:Dtms51:DtmsA + (1  + Dtms + Dtms51 + DtmsA  |Dsu)  ')


% 
Model = fitlme(DATA, 'Drt ~ Dseq + Dseq:Dtms + Dseq:Dtms51  + Dseq:Dtms:DtmsA + Dseq:Dtms51:DtmsA+  Dtms + Dtms51 + DtmsA + Dtms:DtmsA  + Dtms:Dtms51:DtmsA + (1  + Dtms + Dtms51 + DtmsA + Dseq  |Dsu)  ')
%Model = fitlme(DATA, 'Drt ~ Dseq + Dseq:Dtms + Dseq:Dtms51  + Dseq:Dtms:DtmsA + Dseq:Dtms51:DtmsA+  Dtms + Dtms51 + DtmsA + Dtms:DtmsA  + Dtms:Dtms51:DtmsA + (1  + Dseq:Dtms51:DtmsA |Dsu)  ')



% 
Model = fitlme(DATA, 'Drt ~ Dseq + Dseq:Dtms + Dseq:Dtms51  + Dseq:Dtms:DtmsA + Dseq:Dtms51:DtmsA+  Dtms + Dtms51 + DtmsA + Dtms:DtmsA  + Dtms:Dtms51:DtmsA + (1  + Dseq + Dseq:Dtms + Dseq:Dtms51  + Dseq:Dtms:DtmsA + Dseq:Dtms51:DtmsA+  Dtms + Dtms51 + DtmsA + Dtms:DtmsA  + Dtms:Dtms51:DtmsA  |Dsu)  ')

%
Model = fitlme(DATA, 'Drt ~ Dseq + Dseq:Dtms + Dseq:Dtms51  + Dseq:Dtms:DtmsA + Dseq:Dtms51:DtmsA +  Dtms + Dtms51 + DtmsA + Dtms:DtmsA  + Dtms:Dtms51:DtmsA + (1  + Dtms + Dtms51 + DtmsA + Dseq  |Dsu)  ')
 %   AIC           BIC           LogLikelihood    Deviance  
 %   4.8088e+05    4.8112e+05    -2.4041e+05      4.8083e+05

DATAm =DATA;
DATAm.Uno=  DATAm.Dseq==1;
DATA.Uno=  DATA.Dseq==1;
DATAm.Dseq= 1-(1./(DATA.Dseq(:)));


% 
Model = fitlme(DATAm, 'Drt ~ Dseq + Dseq:Dtms + Dseq:Dtms51  + Dseq:Dtms:DtmsA + Dseq:Dtms51:DtmsA +  Dtms + Dtms51 + DtmsA + Dtms:DtmsA  + Dtms:Dtms51:DtmsA + (1  + Dtms + Dtms51 + DtmsA + Dseq  |Dsu)  ')
%   AIC           BIC           LogLikelihood    Deviance   
% 4.8078e+05    4.8101e+05    -2.4036e+05      4.8073e+05


% FINALLLL
Model = fitlme(DATAm, 'Drt ~ Dseq + Dseq:Dtms + Dseq:Dtms51  + Dseq:Dtms:DtmsA + Dseq:Dtms51:DtmsA +  Dtms + Dtms51 + DtmsA +  (1  + Dtms + Dtms51 + DtmsA + Dseq  |Dsu)  ')
%   AIC           BIC           LogLikelihood    Deviance   
% 4.8078e+05    4.8099e+05    -2.4036e+05      4.8073e+05




% 
Model = fitlme(DATAm, 'Drt ~ Dseq + Uno + Dseq:Dtms + Dseq:Dtms51  + Dseq:Dtms:DtmsA + Dseq:Dtms51:DtmsA +  Dtms + Dtms51 + DtmsA +  (1  + Dtms + Dtms51 + DtmsA + Dseq + Uno |Dsu)  ')
%   AIC           BIC           LogLikelihood    Deviance   
%  4.8071e+05    4.8099e+05    -2.4032e+05      4.8065e+05


% 
Model = fitlme(DATAm, 'Drt ~ Dseq + Uno:pError + pError + Dseq:Dtms + Dseq:Dtms51  + Dseq:Dtms:DtmsA + Dseq:Dtms51:DtmsA +  Dtms + Dtms51 + DtmsA +  (1  + Dtms + Dtms51 + DtmsA + Dseq + Uno:pError + pError |Dsu)  ')
%   4.8037e+05    4.8071e+05    -2.4014e+05      4.8029e+05 
%  
% 

% Este parece mas parcimonioso 
Model = fitlme(DATAm, 'Drt ~ Dseq + Uno:pError + Uno + pError + Dseq:Dtms + Dseq:Dtms51  + Dseq:Dtms:DtmsA + Dseq:Dtms51:DtmsA +  Dtms + Dtms51 + DtmsA +  (1  + Dtms + Dtms51 + DtmsA + Dseq + pError  + Uno |Dsu)  ')
%   4.8037e+05    4.8071e+05    -2.4014e+05      4.8029e+05 (x- 1/x)
%   4.8058e+05    4.8094e+05    -2.4025e+05      4.805e+05  (ln(x))
%  4.8059e+05    4.8095e+05    -2.4025e+05      4.8051e+05   exp(-1/x)

%DATAm.Dseq = log (DATA.Dseq);
%DATAm.Dseq= 1-(1./(DATA.Dseq(:)));
%DATAm.Dseq= exp(-1./(DATA.Dseq-1));
Q=0.25;
DATAm.Dseq= 1 -  (1 - Q).^(DATA.Dseq-1);

DATAm.Dseq (DATA.Dseq==1)=0;
DATAm.Dseq (DATA.Dseq==2 | DATA.Dseq==3)=0.25;
DATAm.Dseq ( DATA.Dseq==3)=0.275;
DATAm.Dseq (DATA.Dseq==4| DATA.Dseq==5)=0.3;
DATAm.Dseq (DATA.Dseq==5)=0.4;
DATAm.Dseq (DATA.Dseq==6| DATA.Dseq==7)=0.5;
DATAm.Dseq ( DATA.Dseq==7)=0.75;
% control 1
Model = fitlme(DATAm, 'Drt ~ Dseq + Uno:pError + Uno + pError + pError:Dtms + pError:Dtms51  + pError:Dtms:DtmsA + pError:Dtms51:DtmsA +  Dtms + Dtms51 + DtmsA +  (1  + Dtms + Dtms51 + DtmsA + Dseq + pError  + Uno |Dsu)  ')
% 4.8058e+05    4.8093e+05    -2.4025e+05      4.805e+05  
 
% control 2
Model = fitlme(DATAm, 'Drt ~ Dseq + Uno:pError + Uno + pError + Uno:Dtms + Uno:Dtms51  + Uno:Dtms:DtmsA + Uno:Dtms51:DtmsA +  Dtms + Dtms51 + DtmsA +  (1  + Dtms + Dtms51 + DtmsA + Dseq + pError  + Uno |Dsu)  ')
% 4.806e+05    4.8095e+05    -2.4026e+05      4.8051e+05
 
% control 3
Model = fitlme(DATAm, 'Drt ~ Dseq + Uno:pError + Uno + pError + Uno:pError:Dtms + Uno:pError:Dtms51  + Uno:pError:Dtms:DtmsA + Uno:pError:Dtms51:DtmsA +  Dtms + Dtms51 + DtmsA +  (1  + Dtms + Dtms51 + DtmsA + Dseq + pError  + Uno |Dsu)  ')
% 4.8057e+05    4.8093e+05    -2.4025e+05      4.8049e+05

% control 3bis
Model = fitlme(DATAm, 'Drt ~ Dseq + Uno:pError + Uno + pError + Dseq:Dtms + Dseq:Dtms51  + Dseq:Dtms:DtmsA + Dseq:Dtms51:DtmsA +  Dtms + Dtms51 + DtmsA + Uno:pError:Dtms:DtmsA +  (1  + Dtms + Dtms51 + DtmsA + Dseq + pError  + Uno |Dsu)  ')
% 4.8057e+05    4.8093e+05    -2.4024e+05      4.8048e+05 

% control 3bisbis
Model = fitlme(DATAm, 'Drt ~ Dseq + Uno:pError + Uno + pError + Dseq:Dtms + Dseq:Dtms51  + Dseq:Dtms:DtmsA + Dseq:Dtms51:DtmsA +  Dtms + Dtms51 + DtmsA +  Uno:pError:Dtms + Uno:pError:Dtms51  + Uno:pError:Dtms:DtmsA + Uno:pError:Dtms51:DtmsA +  (1  + Dtms + Dtms51 + DtmsA + Dseq + pError  + Uno |Dsu)  ')
% 4.8057e+05    4.8096e+05    -2.4024e+05      4.8048e+05

% control 3bisbisbis
Model = fitlme(DATAm, 'Drt ~ Dseq + Uno:pError + Uno + pError + Dseq:Dtms + Dseq:Dtms51  + Dseq:Dtms:DtmsA + Dseq:Dtms51:DtmsA +  Dtms + Dtms51 + DtmsA +  Uno:pError:Dtms + Uno:pError:Dtms52  + Uno:pError:Dtms:DtmsA + Uno:pError:Dtms52:DtmsA +  (1  + Dtms + Dtms51 + DtmsA + Dseq + pError  + Uno |Dsu)  ')
% 4.8057e+05    4.8096e+05    -2.4024e+05      4.8048e+05





Model = fitlme(DATAm, 'Drt ~ Dseq + Dseq:Dtms + Dseq:Dtms51  + Dseq:Dtms:DtmsA + Dseq:Dtms51:DtmsA +  Dtms + Dtms51 + DtmsA + Dtms:DtmsA  + Dtms:Dtms51:DtmsA + (1   + Dseq + Dseq:Dtms + Dseq:Dtms51  + Dseq:Dtms:DtmsA + Dseq:Dtms51:DtmsA +  Dtms + Dtms51 + DtmsA + Dtms:DtmsA  + Dtms:Dtms51:DtmsA   |Dsu)  ')
%   AIC           BIC           LogLikelihood    Deviance   
%   4.8078e+05    4.8146e+05    -2.4031e+05      4.8063e+05



DATA_A = DATA(DATAm.DtmsA==1,:);
Model = fitlme(DATA_A, 'Drt ~ Dseq +  Dtms + Dtms51 + Dseq:Dtms + Dseq:Dtms51 +  (1 +Dseq +  Dtms + Dtms51  |Dsu)  ')

DATA_B = DATA(DATAm.DtmsB==1,:);
Model = fitlme(DATA_B, 'Drt ~ Dseq +  Dtms + Dtms51 + Dseq:Dtms + Dseq:Dtms51 +  (1 +Dseq +  Dtms + Dtms51  |Dsu)  ')


%+ Dseq + Dseq:Dtms + Dseq:Dtms51  +  Dtms + Dtms51   + Dtms:Dtms51 
% % second level fixed effect 
% n=0
% clear BE TS
% for s = unique(DATA.Dsu(DATA.DtmsA==0))'
%    n=n+1; 
%   D = DATA(ifcellis(DATA.Dsu,s{1}),:); 
%   R =  fitglm (D,'Drt ~ Dseq + Dseq:Dtms + Dseq:Dtms51  +   Dtms + Dtms51 ') 
%   BE(n) =  table2array(R.Coefficients(5,1))  
%   TS(n) =  table2array(R.Coefficients(5,3))  
% end
% 



%Model = fitlme(DATA, 'Drt ~ Dseq + Dseq*Dtms51 + Dseq*Dtms50 + (1 + Dseq |Dsu)  ')

% esta es!!!!
Model = fitlme(DATA, 'Drt ~     Dtms + Dtms51 + DtmsA + Dtms:DtmsA + DtmsA:Dtms51 + (1 |Dsu)  ')
%98599    98642    -49294           98587   


 
 %%
 
 

%% fix proble in early version of RL.sce

if (sum(RT.est==30)+sum(RT.est==31))==0
    estF =RT.OTHER.F;
    estF(estF==1)=31;estF(estF==0)=30;
    RT.est(3:3:end) = estF(3:3:end); 
    
end

%% fitted a simple Reward lera ning for estract Predicition error and learning signal
%  using JAGS 






% Drt ~ Dseq + Dseq:Dtms + Dseq:Dtms51  + Dseq:Dtms:DtmsA + Dseq:Dtms51:DtmsA +  
%        Dtms + Dtms51 + DtmsA + Dtms:DtmsA  + Dtms:Dtms51:DtmsA 

% Drt ~ Dseq + Uno:pError + Uno + pError + Dseq:Dtms (0.8) + Dseq:Dtms51 (0.079)
% + Dseq:Dtms:DtmsA (0.0173) + Dseq:Dtms51:DtmsA (0) +  Dtms + Dtms51 + DtmsA 


nchains  = 3; % How Many Chains?
nburnin  = 1000; % How Many Burn-in Samples?
nsamples = 1000;  % How Many Recorded Samples?


% Create a single structure that has the data for all observed JAGS nodes
DJ.nT = numel(DATAm.Dseq);
DJ.Dseq = DATAm.Dseq;
DJ.Uno = DATAm.Uno;
DJ.pError = DATAm.pError;
DJ.Dtms = DATAm.Dtms;
DJ.Dtms51 = DATAm.Dtms51;
DJ.DtmsA = DATAm.DtmsA;
DJ.Drt = DATAm.Drt;

su =unique((Dsu(:)));
DJ.idSub=0;
n=0;
warning off
for e =su'
    n=n+1;
   %DJ.Dsu(fun_in_cell(Dsu, [ 'strcmp(@,'  e{1} ') '])) = n;
   ind =  cellfun(@strcmp, Dsu,repmat(e,[ 41245 1]));
   DJ.idSub(ind)=n;
   sum(ind);
end

DJ.Nsubj = n;

% Set initial values for latent variable in each chain
init0=[];
init0(1).mubeta0=0;init0(2).mubeta0=-1;init0(3).mubeta0=1;
init0(1).mubeta1=0;init0(2).mubeta1=1;init0(3).mubeta1=1; 
init0(1).mubeta2=0;init0(2).mubeta2=-1;init0(3).mubeta2=1;
init0(1).mubeta3=0;init0(2).mubeta3=1;init0(3).mubeta3=1; 
init0(1).mubeta4=0;init0(2).mubeta4=-1;init0(3).mubeta4=1;
init0(1).mubeta5=0;init0(2).mubeta5=1;init0(3).mubeta5=1; 
init0(1).mubeta6=0;init0(2).mubeta6=-1;init0(3).mubeta6=1;
init0(1).mubeta7=0;init0(2).mubeta7=1;init0(3).mubeta7=1; 
init0(1).mubeta8=0;init0(2).mubeta8=-1;init0(3).mubeta8=1;
init0(1).mubeta9=0;init0(2).mubeta9=1;init0(3).mubeta9=1; 
init0(1).mubeta10=0;init0(2).mubeta10=-1;init0(3).mubeta10=1;
init0(1).mubeta11=0;init0(2).mubeta11=1;init0(3).mubeta11=1; 
init0(1).mua=0.5;init0(2).mua=2;init0(3).mua=1; 
% fitting
%  U[i] <- V.1[i]  * Pi.1[i]  - V.2[i] * Pi.2[i]
tic
doparallel = 1;
fprintf( 'Running JAGS\n' );
[samples, stats ] = matjags( ...
        DJ, ...
        fullfile('/Volumes/GoogleDrive-112808863907079649330/Mi unidad/Working papers/Martinez_theta/OLDs/DATA', 'model_mixed.txt'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', 10, ...
        'monitorparams', {'mubeta0','mubeta1', 'mubeta2','mubeta3', 'mubeta4','mubeta5', 'mubeta6','mubeta7', 'mubeta8','mubeta9', 'mubeta10','mubeta11'}, ...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0  );
toc

 figure
 subplot(2,2,1)
 histogram(samples.mubeta8(:) )
 title('TMS:Theta:A')
 xlim([-5 22])
 line([ 0 0 ], [ 0 250], 'Color', 'Red','LineWidth',4)
 subplot(2,2,2)
 histogram(samples.mubeta7(:))
  title('TMS:A')
 line([ 0 0 ], [ 0 250], 'Color', 'Red','LineWidth',4)
  subplot(2,2,3)
 histogram(samples.mubeta6(:))
  title('TMS:Theta')
 line([ 0 0 ], [ 0 250], 'Color', 'Red','LineWidth',4) 
  subplot(2,2,4)
 histogram(samples.mubeta5(:))
  title('TMS')
 line([ 0 0 ], [ 0 250], 'Color', 'Red','LineWidth',4) 

 
 
 
 %LUGAR = '/Users/pablobilleke/Desktop/_______________PASO____________________';
LUGAR = '/Volumes/GoogleDrive-112808863907079649330/Mi unidad/Working papers/Martinez_theta/OLDs/DATA';

 tic
doparallel = 1;
fprintf( 'Running JAGS\n' );
[samplesR, statsR ] = matjags( ...
        DJ, ...
        fullfile(LUGAR, 'model_mixed_robust.txt'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', 10, ...
        'monitorparams', {'mubeta0','mubeta1', 'mubeta2','mubeta3', 'mubeta4','mubeta5', 'mubeta6','mubeta7', 'mubeta8','mubeta9', 'mubeta10','mubeta11'}, ...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0  );
toc


 figure
 subplot(2,2,1)
 histogram(samplesR.mubeta8(:) )
 title('TMS:Theta:A')
 xlim([-5 22])
 line([ 0 0 ], [ 0 250], 'Color', 'Red','LineWidth',4)
 subplot(2,2,2)
 histogram(samplesR.mubeta7(:))
  title('TMS:A')
 line([ 0 0 ], [ 0 250], 'Color', 'Red','LineWidth',4)
  subplot(2,2,3)
 histogram(samplesR.mubeta6(:))
  title('TMS:Theta')
 line([ 0 0 ], [ 0 250], 'Color', 'Red','LineWidth',4) 
  subplot(2,2,4)
 histogram(samplesR.mubeta5(:))
  title('TMS')
 line([ 0 0 ], [ 0 250], 'Color', 'Red','LineWidth',4) 

 %%
 
  %LUGAR = '/Users/pablobilleke/Desktop/_______________PASO____________________';
LUGAR = '/Volumes/GoogleDrive-112808863907079649330/Mi unidad/Working papers/Martinez_theta/OLDs/DATA';

DJ.Dseq = DATA.Dseq;

 tic
doparallel = 1;
fprintf( 'Running JAGS\n' );
[samplesMA, statsMA ] = matjags( ...
        DJ, ...
        fullfile(LUGAR, 'model_mixed_alpha.txt'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', 10, ...
        'monitorparams', {'mubeta0','mubeta1', 'mubeta2','mubeta3', 'mubeta4','mubeta5', 'mubeta6','mubeta7', 'mubeta8','mubeta9', 'mubeta10','mubeta11','mua'}, ...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0  );
toc

  figure
 subplot(2,2,1)
 histogram(samplesMA.mubeta8(:) )
 title('TMS:Theta:A')
 xlim([-5 22])
 line([ 0 0 ], [ 0 250], 'Color', 'Red','LineWidth',4)
 subplot(2,2,2)
 histogram(samplesMA.mubeta7(:))
  title('TMS:A')
 line([ 0 0 ], [ 0 250], 'Color', 'Red','LineWidth',4)
  subplot(2,2,3)
 histogram(samplesMA.mubeta6(:))
  title('TMS:Theta')
 line([ 0 0 ], [ 0 250], 'Color', 'Red','LineWidth',4) 
  subplot(2,2,4)
 histogram(samplesMA.mubeta5(:))
  title('TMS')
 line([ 0 0 ], [ 0 250], 'Color', 'Red','LineWidth',4) 

 
%%

nchains  = 3; % How Many Chains?
nburnin  = 1000; % How Many Burn-in Samples?
nsamples = 1000;  % How Many Recorded Samples?


% Create a single structure that has the data for all observed JAGS nodes
DATA.Uno=  DATA.Dseq==1;
DJ.nT = numel(DATA.Dseq);
DJ.Dseq = DATA.Dseq;
DJ.Uno = DATA.Uno;
DJ.pError = DATA.pError;
DJ.Dtms = DATA.Dtms;
DJ.Dtms51 = DATA.Dtms51;
DJ.DtmsA = DATA.DtmsA;
DJ.Drt = DATA.Drt;

su =unique((Dsu(:)));
DJ.idSub=0;
n=0;
warning off
for e =su'
    n=n+1;
   %DJ.Dsu(fun_in_cell(Dsu, [ 'strcmp(@,'  e{1} ') '])) = n;
   ind =  cellfun(@strcmp, Dsu,repmat(e,[ 41245 1]));
   DJ.idSub(ind)=n;
   sum(ind);
end

DJ.Nsubj = n;


% Set initial values for latent variable in each chain
init0=[];
init0(1).mubeta0=0;init0(2).mubeta0=-1;init0(3).mubeta0=1;
init0(1).mubeta1=0;init0(2).mubeta1=1;init0(3).mubeta1=1; 
init0(1).mubeta2=0;init0(2).mubeta2=-1;init0(3).mubeta2=1;
init0(1).mubeta3=0;init0(2).mubeta3=1;init0(3).mubeta3=1; 
init0(1).mubeta4=0;init0(2).mubeta4=-1;init0(3).mubeta4=1;

init0(1).muasham =0.3;init0(2).muasham =0.5;init0(3).muasham =0.9; 
init0(1).muathA =0.3;init0(2).muathA =0.5;init0(3).muathA =0.9; 
init0(1).muathB =0.3;init0(2).muathB =0.5;init0(3).muathB =0.9; 
init0(1).muanthA =0.3;init0(2).muanthA =0.5;init0(3).muanthA =0.9; 
init0(1).muanthB =0.3;init0(2).muanthB =0.5;init0(3).muanthB =0.9; 

init0(1).mubeta5=0;init0(2).mubeta5=-1;init0(3).mubeta5=1;
init0(1).mubeta6=0;init0(2).mubeta6=1;init0(3).mubeta6=1; 
init0(1).mubeta7=0;init0(2).mubeta4=-1;init0(3).mubeta7=1;

%LUGAR = '/Users/pablobilleke/Desktop/_______________PASO____________________';
LUGAR = '/Volumes/GoogleDrive-112808863907079649330/Mi unidad/Working papers/Martinez_theta/OLDs/DATA';
tic
doparallel = 0;
fprintf( 'Running JAGS\n' );
[samplesC, statsC ] = matjags( ...
        DJ, ...
        fullfile(LUGAR, 'model_cognitive_proport_.txt'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', 10, ...
        'monitorparams', {'mubeta0','mubeta1', 'mubeta2','mubeta3', 'mubeta4','muasham', 'muathA','muathB', 'muanthA','muanthB','mubeta5','mubeta6', 'mubeta7'}, ...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0  );
toc




%% Accuaracy


%%

%
Error  = (TRT.rt>10 & TRT.est==21) | (TRT.rt<0 & TRT.est==10);
NGO = TRT.est==21;
pError = [ 0  Error(1:end-1) ] ;
Drt = TRT.rt;%( TRT.rt>10 & TRT.est==10 )';
pRT = [ -99 Drt(1:end-1) ] ;
paso = find(TRT.OTHER.seq==0);
TRT.OTHER.seq(paso) = TRT.OTHER.seq(paso-1);
Dseq =TRT.OTHER.seq;
DtmsB = ifcellis(  TRT.OTHER.TMS_SITE, 'B');
DtmsA = ifcellis(  TRT.OTHER.TMS_SITE, 'A');
Dtms = TRT.OTHER.TMS';
Dsu = TRT.OTHER.Sub';

nAC = double(Error(NGO))';
pError =pError(NGO)'; 
Dsu = Dsu(NGO);

pRT =pRT(NGO)';
pRT(pRT<10)=0;
Dseq = Dseq(NGO)';

Dtms = Dtms(NGO);
DtmsB = DtmsB(NGO)';
DtmsA = DtmsA(NGO)';

Dtms51 = Dtms==51;
Dtms52 = Dtms==52;
Dtms50 = Dtms==50;
Dtms=Dtms51+Dtms52;


DATA = table(pRT,Dseq,Dsu,Dtms51,Dtms52,Dtms50,Dtms, DtmsB, DtmsA, nAC, pError);

Model = fitlme(DATA, 'nAC ~ Dseq +  pError + pRT  +  (1  + Dseq +  pError + pRT  |Dsu)  ');

Model = fitlme(DATA, 'nAC ~ Dseq +  pError + pRT  +  Dtms + Dtms51 + DtmsA + Dseq:Dtms  + Dseq:DtmsA + Dseq:Dtms51 + Dseq:Dtms51:DtmsA+  (1  + Dseq +  pError + pRT +  Dtms + Dtms51 + DtmsA |Dsu)  ');

nchains  = 3; % How Many Chains?
nburnin  = 1000; % How Many Burn-in Samples?
nsamples = 1000;  % How Many Recorded Samples?


% Create a single structure that has the data for all observed JAGS nodes
DJ.nT = numel(DATA.Dseq);
DJ.Dseq = DATA.Dseq;
DJ.pError = DATA.pError;
DJ.Dtms = DATA.Dtms;
DJ.Dtms51 = DATA.Dtms51;
DJ.DtmsA = DATA.DtmsA;
DJ.nAC = DATA.nAC;
DJ.pRT = DATA.pRT;

su =unique((Dsu(:)));
DJ.idSub=0;
n=0;
warning off
for e =su'
    n=n+1;
   %DJ.Dsu(fun_in_cell(Dsu, [ 'strcmp(@,'  e{1} ') '])) = n;
   ind =  cellfun(@strcmp, Dsu,repmat(e,[ 10171 1]));
   DJ.idSub(ind)=n;
   sum(ind);
end

DJ.Nsubj = n;


% Set initial values for latent variable in each chain
init0=[];
init0(1).mubeta0=0;init0(2).mubeta0=-.1;init0(3).mubeta0=.1;
init0(1).mubeta1=0;init0(2).mubeta1=.1;init0(3).mubeta1=.1; 
init0(1).mubeta2=0;init0(2).mubeta2=-.1;init0(3).mubeta2=.1;
init0(1).mubeta3=0;init0(2).mubeta3=.1;init0(3).mubeta3=.1; 
init0(1).mubeta4=0;init0(2).mubeta4=-1;init0(3).mubeta4=1;
init0(1).mubeta5=0;init0(2).mubeta5=1;init0(3).mubeta5=1; 
init0(1).mubeta6=0;init0(2).mubeta6=-1;init0(3).mubeta6=1;
init0(1).mubeta7=0;init0(2).mubeta7=1;init0(3).mubeta7=1; 
init0(1).mubeta8=0;init0(2).mubeta8=-1;init0(3).mubeta8=1;
init0(1).mubeta9=0;init0(2).mubeta9=1;init0(3).mubeta9=1; 
init0(1).mubeta10=0;init0(2).mubeta10=-1;init0(3).mubeta10=1;
%init0(1).mubeta4=0;init0(2).mubeta4=-1;init0(3).mubeta4=1;

%LUGAR = '/Users/pablobilleke/Desktop/_______________PASO____________________';
LUGAR = '/Volumes/GoogleDrive-112808863907079649330/Mi unidad/Working papers/Martinez_theta/OLDs/DATA';
tic
doparallel = 1;
fprintf( 'Running JAGS\n' );
[samplesAC, statsAC ] = matjags( ...
        DJ, ...
        fullfile(LUGAR, 'model_mixed_AC.txt'), ...
        init0, ...
        'doparallel' , doparallel, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', 10, ...
        'monitorparams', {'mubeta0','mubeta1', 'mubeta2','mubeta3','mubeta4','mubeta5', 'mubeta6','mubeta7','mubeta8','mubeta9', 'mubeta10'}, ...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0  );
toc



  figure
 subplot(2,2,1)
 histogram(samplesAC.mubeta0(:) )
 title('int')
 xlim([-5 22])
 line([ 0 0 ], [ 0 250], 'Color', 'Red','LineWidth',4)
 subplot(2,2,2)
 histogram(samplesAC.mubeta1(:))
  title('Seq')
 line([ 0 0 ], [ 0 250], 'Color', 'Red','LineWidth',4)
  subplot(2,2,3)
 histogram(samplesAC.mubeta2(:))
  title('pError')
 line([ 0 0 ], [ 0 250], 'Color', 'Red','LineWidth',4) 
  subplot(2,2,4)
 histogram(samplesAC.mubeta3(:))
  title('pRT')
 line([ 0 0 ], [ 0 250], 'Color', 'Red','LineWidth',4) 
