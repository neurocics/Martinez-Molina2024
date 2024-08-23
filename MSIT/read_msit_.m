clear


PATH = '/Volumes/neuroCICS01/sujetos_tesis_pcp/'


SUB= {'ACRC_11051992'
'AICR_06072000'
'AVOL_02092003'
'B_FE_24011989'
'BAHA_27122000'
'BCCN_26051999'
'BCGG_21031989'
'BMDP_26012004'
'BRSG_12101999'
'CAGG_14031988'
'CAHG_27061988'
'CARA_10081990'
'CBPR_29121988'
'CDPR_25031986'
'CEAS_23071992'
'CFPS_05081994'
'CPAS_09102000'
'CPVA_31031992'
'DGAH_11092000'
'FDCV_15092001'
'GACM_21011991'
'IAGR_02091993'
'IAPR_03051989'
'ICMP_09051991'
'JEYO_29011989'
'JLZB_01031991'
'JMBV_03071992'
'JMPP_04061996'
'JPRC_29061986'
'KECT_07071995'
'KGUR_28072004'
'LFSC_05031993'
'LPLC_08111999'
'LTGS_17051996'
'MAAO_15081988'
'MABS_24121998'
'MEFR_14061991'
'MINR_21112000'
'MJCO_15032000'
'MJVT_03102002'
'MJZM_09031993'
'MMCP_13011990'
'MQAB_30011998'
'NELS_07081998'
'NJRA_24051999'
'NPLM_05071997'
'PACP_05022004'
'PAYO_17041984'
'PDLR_27011999'
'PXGV_07041997'
'RAGG_26081990'
'RJNV_07091989'
'SACC_19121997'
'SNCS_22121989'
'TPPS_14091998'
'YDCL_11021994'}'

for  s = 1:numel(SUB)

%cd(([ PATH  '/' SUB{s}  '/']))

 clear RT

if 1 % calcualr las correctas
file = ls_lan([ PATH  '/' SUB{s}  '/LOG/*fMRI*MSIT_pb.log'])
cfg = [];
cfg.filename =  file{1}
cfg.type = 'presentation' 
cfg.est = [11 12]
cfg.resp = [1 2 3]
cfg.rw = [10000]       %   (ventada de rerspuestas, en cfg.unit)(for TR in MS!!!)
RT = rt_read(cfg)
numel(RT.est)
file = ls([ PATH  '/'  SUB{s}  '/LOG/*fMRI*all_Respuestas_*.txt']);
file(isspace(file))=[];
msit = import_msit(file);

if  strcmp(SUB{s},'CCAC_12071987')

    msit = msit(1:275,:)

end
acu = (((msit.l1==msit.r)+(msit.l2==msit.r)+(msit.l3==msit.r))==1);
RT.OTHER.acu = acu(:)'
RT.OTHER.seq = msit.N(:)'
RT.OTHER.subN = ones(size(RT.OTHER.acu))*s
RT.OTHER.subID = cell(size(RT.OTHER.acu))
RT.OTHER.subID(:)=SUB(s)
save ([ PATH  '/'   SUB{s}  '/LOG/RTmsit' ], 'RT')

end


if s==1
    RTall = RT;
else

RTall = rt_merge(RTall,RT);

end
end


%%

COR=[]
COR.RT=RTall

COR2tableR(COR)



