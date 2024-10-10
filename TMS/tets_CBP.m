
clear
close all
load('TMS_model_22_23s.mat')
%load('model_grupal_exp.mat')

%figure_lan
%topoplot_lan([],GLAN.chanlocs,'style','blank','electrodes','labelpoint');  

show_elec(GLAN)

IZall = [1 51 37 11 45 3 31 53 25 39 21 13 47 5 33];
DERall = IZall +1;
%
t1 = find_approx(GLAN.timefreq.time,-0.5) ;
t2 = find_approx(GLAN.timefreq.time,1) ;


clear DATA
DATA{1} = GLAN.timefreq.subdata{2,5}(:,:,t1:t2,:); 
%DATA{2} = GLAN.timefreq.subdata{2,2}(:,:,t1:t2,:); 

%DATA{1} = (DATA{1} + DATA{2})/2;
%DATA = DATA(1)

[pval , stats] = lan_nonparametric(DATA);


[pvalc cluster] = lan_CBPt(DATA,GLAN.chanlocs(1).electrodemat,0.05,1000,0.5) ;

unique(pvalc)
%%

t1 = find_approx(GLAN.timefreq.time,-0.5) ;
t2 = find_approx(GLAN.timefreq.time,1) ;
f1 = find_approx(GLAN.timefreq.freq,1) ;
f2 = find_approx(GLAN.timefreq.freq,30) ;



ELE=[IZall DERall];
p_plot = pvalc<0.05 & pval<0.05;
p_plot = sum(p_plot(:,ELE,:),2);



w_plot = nanmean(DATA{1},4);
w_plot(pvalc>0.05 & pval>0.05)=NaN;
w_plot = nansum(w_plot(:,ELE,:),2);
w_plot(abs(w_plot)==0 | p_plot<0)=0;
w_plot_h = w_plot ;
w_plot_h(p_plot<1)=NaN;
w_plot_h = reduce_cluster(squeeze(w_plot_h), 50,8);


if 0
 figure_lan
 pcolor2(GLAN.timefreq.time(t1:t2), GLAN.timefreq.freq, squeeze(p_plot )), shading flat
end

if 1
 figure_lan
 H = pcolor2(GLAN.timefreq.time(t1:t2), GLAN.timefreq.freq(f1:f2), squeeze(w_plot )), shading flat
 alpha(H,0.5)  
 hold on 
 pcolor2(GLAN.timefreq.time(t1:t2), GLAN.timefreq.freq(f1:f2), squeeze(w_plot_h )), shading flat
end

%save CPT_seq_TMS_nTHh_B pval* cl* DATA
%save CPT_seq_AB pval* cl*
% save CPT_seq_TMS_THh_A pval* cl* DATA


%%
figure_lan
t1t= find_approx(GLAN.timefreq.time,0.25) ;
t2t = find_approx(GLAN.timefreq.time,0.5) ;
f1t=find_approx(GLAN.timefreq.freq,4.5) ;
f2t=find_approx(GLAN.timefreq.freq,6.5) ;

tdata = squeeze(nanmean(nanmean(nanmean(GLAN.timefreq.subdata{1,2}(f1t:f2t,:,t1t:t2t,:),4),3),1));

tdata = tdata + squeeze(nanmean(nanmean(nanmean(GLAN.timefreq.subdata{2,2}(f1t:f2t,:,t1t:t2t,:),4),3),1));

topoplot_lan(tdata/2,GLAN.chanlocs, 'style','map' ,'electrodes','off','whitebk', 'on' )


