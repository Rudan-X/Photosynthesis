clear model

model=struct();
model.S=zeros(0,0);
model.rxns=cell(0);
model.mets=cell(0);

model.lb=zeros(0);
model.ub=ones(0)*1000;
model.rxnNames=cell(0);
model.b=zeros(0);
model.c=zeros(0);

cat=cell(0); % record which enzyme catalyze which reaction


% Light absorption by the antenna
model=addReaction(model,'import_photon',' -> photon');
model=addReaction(model,'irra_m2','chl2_m + 6 photon -> chl2_me');
model=addReaction(model,'irra_peri2','chl2_p + 6 photon -> chl2_pe');
model=addReaction(model,'irra_core2','chl2_c + 6 photon -> chl2_ce');
model=addReaction(model,'irra_m1','chl1_m + 6 photon -> chl1_me');
model=addReaction(model,'irra_peri1','chl1_p + 6 photon -> chl1_pe');
model=addReaction(model,'irra_core1','chl1_c + 6 photon -> chl1_ce');


% Electron transport through antenna
model=addReaction(model,'m2e_c2','chl2_me + chl2_c <=> chl2_ce + chl2_m');
model=addReaction(model,'p2e_c2','chl2_pe + chl2_c <=> chl2_ce + chl2_p');
model=addReaction(model,'m1e_c1','chl1_me + chl1_c <=> chl1_ce + chl1_m');
model=addReaction(model,'p1e_c1','chl1_pe + chl1_c <=> chl1_ce + chl1_p');

% PSII
model=addReaction(model,'c2e_680','chl2_ce + P680 -> P680e + chl2_c');
model=addReaction(model,'p680e_680p','P680e + pheo -> pheon + P680p');
model=addReaction(model,'pheo_qa','pheon + QA -> pheo + QAn');
model=addReaction(model,'qa_qb','QAn + QB -> QA + QBn');
model=addReaction(model,'qa_qbn','QAn + QBn -> QA + QBn2');
model=addReaction(model,'pq_pqh2','QBn2 + PQ + 2 H[s] -> PQH2 + QB');


% Oxygen evolving complex giving electron to RC, all irreversible
model=addReaction(model,'s0_s0p','P680p + s0 -> P680 + s0p');
model=addReaction(model,'s0p_s1','s0p -> s1');

model=addReaction(model,'s1_s1p','P680p + s1 -> P680 + s1p');
model=addReaction(model,'s1p_s2','s1p -> s2');

model=addReaction(model,'s2_s2p','P680p + s2 -> P680 + s2p');
model=addReaction(model,'s2p_s3','s2p -> s3');

model=addReaction(model,'s3_s3p','P680p + s3 -> P680 + s3p');
model=addReaction(model,'s3p_s4','s3p -> s4');

model=addReaction(model,'s4_s4p','2 H2O[l] + s4 -> O2[l] + 4 H[l] + s4p');
model=addReaction(model,'s4p_s0','s4p -> s0');

%% Cytochrome
model=addReaction(model,'pqh2_ispox','PQH2 + ISPox -> ISPoxPQH2');
model=addReaction(model,'ispoxp','ISPoxPQH2 -> PQsemi + ISPHred');
model=addReaction(model,'isphr','ISPHred + cytf -> cytfn + ISPHp');
model=addReaction(model,'isphp','ISPHp <=> ISPox + H[l]');
model=addReaction(model,'cytf_pc','PC + cytfn <=> cytf + PCn');
model=addReaction(model,'pqse_pq','PQsemi + cytbl -> PQ + cytbln + H[l]'); 
model=addReaction(model,'cytl_h','cytbln + cytbh -> cytbl + cytbhn');
model=addReaction(model,'cyt_pq','cytbhn + PQ -> cytbh + PQn');
model=addReaction(model,'cyt_pq2','cytbhn + PQn -> cytbh + PQn2');
model=addReaction(model,'pqproto','PQn2 + 2 H[s] -> PQH2');


% PSI
% Electron transport through antenna
model=addReaction(model,'c1e_700','chl1_ce + P700 -> P700e + chl1_c');
model=addReaction(model,'p700e_700p','P700e + A0 -> A0n + P700p');
model=addReaction(model,'p700p_700','P700p + PCn -> P700 + PC');
model=addReaction(model,'a_fd','A0n + Fd -> Fdn + A0');

% Ion movement
%model=addReaction(model,'mgsl','mg2p[s] <=> mg2p[l]'); 
model=addReaction(model,'ksl','k[s] <=> k[l]'); 


% ATP activase activation and ATP NADPH synthesis
model=addReaction(model,'atps_ac','atps + Fdn <=> atpsac + Fd'); 
model=addReaction(model,'atp','ADP[s] + Pi[s] + 4 H[l] <=> ATP[s] + 4 H[s] + H2O[s]'); cat(length(model.rxns))={'atpsac'};
model=addReaction(model,'fnr','2 Fdn + NADPp[s] + H[s] -> 2 Fd + NADPH[s]'); cat(length(model.rxns))={'fnr'}; 

% Thio activation
model=addReaction(model,'ftr','2 Fdn + thio <=> 2 Fd + thion + 2 H[s]'); cat(length(model.rxns))={'ftr'};

% Enzyme activation
model=addReaction(model,'rca_ac','Rcai + ATP[s] -> Rcaa + ADP[s]');
model=addReaction(model,'RA1','ER + Rcaa -> ERRcaa');
model=addReaction(model,'RA2','ERRcaa -> RuBP + ERcaa');
model=addReaction(model,'RA3','ERcaa -> E + Rcai');

model=addReaction(model,'RA1b','ECR + Rcaa -> ECRRcaa');
model=addReaction(model,'RA2b','ECRRcaa -> RuBP + ECRcaa');
model=addReaction(model,'RA3b','ECRcaa -> EC + Rcai');

model=addReaction(model,'EC','E + CO2[s] <=> EC'); 
model=addReaction(model,'ECM','EC + mg2p[s] <=> ECM'); 
model=addReaction(model,'ER','E + RuBP -> ER');
model=addReaction(model,'ECR','EC + RuBP -> ECR'); 
model=addReaction(model,'ECMR','ECM + RuBP <=> ECMR');

model=addReaction(model,'carb1','ECMR + CO2[s] -> ECMRC');
model=addReaction(model,'carb2','ECMRC -> RPGA2H');
model=addReaction(model,'carb3','RPGA2H -> RPGA2 + H[s]');
model=addReaction(model,'carb4','RPGA2 -> RPGA + PGA[s]');
model=addReaction(model,'carb5','RPGA -> ECM + PGA[s]');

% RuBP can be produced during nucleotide metabolism 

%%
% Calvin Cycle
model=addReaction(model,'pgaks','ATP[s] + PGA[s] <=> ADP[s] + DpGA[s]'); cat(length(model.rxns))={'pgak_s'}; % PGA Kinase 2.7.2.3	PGA+ATP -> ADP + DPGA
model=addReaction(model,'fbpasen','FBP[s] + H2O[s] -> F6P[s] + Pi[s]'); cat(length(model.rxns))={'fbpasen'}; % FBPasen(fbpase) 3.1.3.11	FBP->F6P+OP
model=addReaction(model,'fbpasenb','F6P[s] + Pi[s] -> FBP[s] + H2O[s]'); cat(length(model.rxns))={'fbpase'};
model=addReaction(model,'gapdhs','DpGA[s] + NADPH[s] <=> NADPp[s] + GAP[s] + Pi[s]'); cat(length(model.rxns))={'gapdh_s'}; % GAP dehydragenase (GAPDH) 1.2.1.13	DPGA+NADPH <->GAP + OP+NADP 
model=addReaction(model,'tpis','DHAP[s] <=> GAP[s]'); cat(length(model.rxns))={'tpi_s'}; % Triose phosphate isomerase 5.3.1.1	DHAP <->GAP
model=addReaction(model,'fbpaldos','GAP[s] + DHAP[s] <=> FBP[s]'); cat(length(model.rxns))={'fbpaldo_s'}; % Aldolase	4.1.2.13 GAP+DHAP <->FBP


model=addReaction(model,'transk','F6P[s] + GAP[s] <=> E4P + X5P'); cat(length(model.rxns))={'transk'}; % Transketolase 2.2.1.1 F6P+GAP<->E4P+Xu5P
model=addReaction(model,'fbpaldos2','E4P + DHAP[s] <=> SBP'); cat(length(model.rxns))={'fbpaldo_s'}; % Aldolase(aldoa) 4.1.2.13	E4P+DHAP<->SBP
model=addReaction(model,'sbpasen','SBP + H2O[s] -> S7P + Pi[s]'); cat(length(model.rxns))={'sbpasen'}; % SBPasen(sbpase) 3.1.3.37	SBP->S7P+OP
model=addReaction(model,'sbpasenb','S7P + Pi[s] -> SBP + H2O[s]'); cat(length(model.rxns))={'sbpase'}; 


model=addReaction(model,'transk2','S7P + GAP[s] <=> R5P + X5P'); cat(length(model.rxns))={'transk'}; % Transketolase(tk) 2.2.1.1	S7P+GAP<->Ri5P+X5P
model=addReaction(model,'rpi','R5P <=> RU5P'); cat(length(model.rxns))={'rpi'};  %% % Pentosephosphate isomerase(rpi) 5.3.1.6 R5P<-->Ru5P
model=addReaction(model,'rpe','X5P <=> RU5P'); cat(length(model.rxns))={'rpe'}; % Pentosephosphate epimerase(rpe) 5.1.3.1 Xu5P<-->Ru5P
model=addReaction(model,'prkn','RU5P + ATP[s] -> RuBP + ADP[s]'); cat(length(model.rxns))={'prkn'}; % PRKn	2.7.1.19 Ru5P+ATP ->RuBP+ADP
model=addReaction(model,'prknb','RuBP + ADP[s] -> RU5P + ATP[s]'); cat(length(model.rxns))={'prk'}; 

% synthesis of starch in chloroplast
model=addReaction(model,'gpis','F6P[s] <=> G6P[s]'); cat(length(model.rxns))={'gpi_s'}; % Hexose phosphate isomerase(gpi) 5.3.1.9	F6P<->G6P
model=addReaction(model,'pgms','G6P[s] <=> G1P[s]'); cat(length(model.rxns))={'pgm_s'}; % Phosphoglucomutase(pgm) 5.4.2.2 G6P<->G1P
model=addReaction(model,'adpgpy','G1P[s] + ATP[s] + H[s] -> ADPG + Pi[s]'); cat(length(model.rxns))={'adpgpy'}; % ADP-glucose pyrophosphorylase 2.7.7.27: G1P+ATP ->ADPG+ADP
model=addReaction(model,'ss','ADPG + glucan -> STA + ADP[s]'); cat(length(model.rxns))={'starsyn'}; % Starch synthase (2.4.1.21), 1,4-alpha- branching enzyme (2.4.1.18) ADPG+Glucan -> STA 
model=addReaction(model,'agb','STA + Pi[s] -> glucan + G1P[s]'); cat(length(model.rxns))={'agb'}; %1,4-alpha-glucan branching enzyme 2.4.1.1  ADPG+Gn<->G(n+1): STA + Pi -> Glucan + G1P
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1081201/
% added import of ADPG from cytosol to chloroplast
model=addReaction(model,'adpgimp',' -> ADPG'); 


% exchange of DHAP, GAP, PGA and DPGA between chloroplast and stroma
model=addReaction(model,'vdhap_in','DHAP[s] + Pi[c] <=> DHAP[c] + Pi[s]'); 
model=addReaction(model,'vgap_in','GAP[s] + Pi[c] <=> GAP[c] + Pi[s]'); 
model=addReaction(model,'vpga_in','PGA[s] + Pi[c] <=> PGA[c] + Pi[s]'); 
model=addReaction(model,'vdpga_in','DpGA[s] + Pi[c] <=> DpGA[c] + Pi[s]'); 


% Phosphorylation
model=addReaction(model,'resp1','ECMR + O2[s] -> ECMRO2');
model=addReaction(model,'resp2','ECMRO2 -> PGAPGLRubiH2');
model=addReaction(model,'resp3','PGAPGLRubiH2 -> PGAPGLRubi + 2 H[s]');
model=addReaction(model,'resp4','PGAPGLRubi -> PGAPGL + ECM');
model=addReaction(model,'resp5','PGAPGL -> PGA[s] + PGL');

% Photorespiration
model=addReaction(model,'gca','2 PGL + H2O[s] -> GCA[s]');
model=addReaction(model,'gcaex','GCA[s] <=> GCA[c]');
model=addReaction(model,'goa','GCA[c] + O2[c] -> GOA + H2O2[c]');
model=addReaction(model,'h2o2trans','2 H2O2[c] -> O2[c] + 2 H2O[c]');
model=addReaction(model,'gly','GOA <=> GLY');
model=addReaction(model,'ser','2 GLY + NADp -> SER + CO2[c] + NADH');
model=addReaction(model,'hpr','SER + GOA <=> HPR + GLY');
model=addReaction(model,'gceac','HPR + NADH <=> GCEA[c] + NADp');
model=addReaction(model,'gceaex','GCEA[c] <=> GCEA[s]');
model=addReaction(model,'pgareg','GCEA[s] + ATP[s] <=> PGA[s] + ADP[s]');
model=addReaction(model,'pgasource',' <=> PGA[s]');

% Synthesis of glucose im cytosol 
model=addReaction(model,'pgakc','PGA[c] + ATP[c] <=> ADP[c] + DpGA[c]'); cat(length(model.rxns))={'pgak_c'}; % PGA Kinase 2.7.2.3	PGA+ATP -> ADP + DPGA
model=addReaction(model,'gapdhc','DpGA[c] + H[c] <=> GAP[c] + Pi[c]'); cat(length(model.rxns))={'gapdh_c'}; % GAP dehydragenase (GAPDH) 1.2.1.13	DPGA+NADPH <->GAP + OP+NADP 
model=addReaction(model,'tpic','DHAP[c] <=> GAP[c]'); cat(length(model.rxns))={'tpi_c'}; % Triose phosphate isomerase 5.3.1.1	DHAP <->GAP
model=addReaction(model,'fbpaldoc','GAP[c] + DHAP[c] <=> FBP[c]'); cat(length(model.rxns))={'fbpaldo_c'}; % Aldolase	4.1.2.13 GAP+DHAP <->FBP
model=addReaction(model,'fbpasec','FBP[c] + H2O[c] <=> F6P[c] + Pi[c]'); cat(length(model.rxns))={'fbpase_c'}; % FBPase(fbpase) 3.1.3.11	FBP<->F6P+OP

model=addReaction(model,'gpic','F6P[c] <=> G6P[c]'); cat(length(model.rxns))={'gpi_c'}; % Hexose phosphate isomerase(gpi) 5.3.1.9	F6P<->G6P
model=addReaction(model,'pgmc','G6P[c] <=> G1P[c]'); cat(length(model.rxns))={'pgm_c'}; % Phosphoglucomutase(pgm) 5.4.2.2 G6P<->G1P
model=addReaction(model,'ugpase','G1P[c] + UTP + H[c] <=> 2 Pi[c] + UDPG'); cat(length(model.rxns))={'ugpase'}; % UDP glucose pyrohpsphorylase(ugpase) 2.7.7.9: G1P+UTP -->2Pi +UDPG 
model=addReaction(model,'sps','UDPG + F6P[c] <=> SUCP + UDP'); cat(length(model.rxns))={'sps'};% sucrose phosphate synthase(sps) (2.4.1.14) UDPG+F6P-->SUCP + UDP + H + Pi
model=addReaction(model,'spase','SUCP <=> SUC + Pi[c]'); cat(length(model.rxns))={'spase'}; % sucrose phosphatase phosphatase(spp) (3.1.3.24): SUCP--Pi + SUC
model=addReaction(model,'ndpk','ATP[c] + UDP <=> UTP + ADP[c]'); cat(length(model.rxns))={'ndpk'}; % nucleoside-diphosphate kinase(ndpk) 2.7.4.6: ATP+UDP --UTP + ADP
model=addReaction(model,'f26bpase','F26BP <=> Pi[c] + F6P[c]'); cat(length(model.rxns))={'f26bpase'}; % fructose-2,6-bisphosphate 2-phosphatase (f26bp) 3.1.3.46: F26BP--F6P + Pi
model=addReaction(model,'pfk2','F6P[c] + ATP[c] <=> ADP[c] + F26BP'); cat(length(model.rxns))={'pfk2'}; % 6-phosphofructo-2-kinase(pfk2) 2.7.1.105: F6P + ATP --ADP + F26BP



% Regulation

% high light
model=addReaction(model,'fbpase','fbpase + thion <=> fbpasen + thio');
model=addReaction(model,'sbpase','sbpase + thion <=> sbpasen + thio');
model=addReaction(model,'prk','prk + thion <=> prkn + thio');
% model=addReaction(model,'cp12','CP + thion <=> CPn + thio');
% model=addReaction(model,'cgp','CP + PRK + GAPDH[s] <=> CGP');


% Darkness
model=addReaction(model,'fbpase_in','fbpasen + Trxl2 -> fbpase + Trxl2n');
model=addReaction(model,'sbpase_in','sbpasen + Trxl2 -> sbpase + Trxl2n');
model=addReaction(model,'prk_in','prkn + Trxl2 -> prk + Trxl2n');
% model=addReaction(model,'cp12_in','CPn + Trxl2 -> CP + Trxl2n');
% model=addReaction(model,'Rcaa_in','Rcaa + Trxl2 -> Rcai + Trxl2n');

model=addReaction(model,'on2','Fdn + O2[s] -> Fd + On2');
model=addReaction(model,'h2o2','2 On2 + 2 H[s] -> H2O2[s] + O2[s]');
model=addReaction(model,'cysprx','H2O2[s] + cysprxn -> H2O[s] + cysprx');
model=addReaction(model,'trxl2','Trxl2n + cysprx -> Trxl2 + cysprxn');


% Photochemical quenching
model=addReaction(model,'psbs_act','psbs + 6 H[l] <=> psbsh6');
model=addReaction(model,'vde_act','VDE + 3 H[l] <=> vdeh3');
model=addReaction(model,'va','Vx -> Ax'); cat(length(model.rxns))={'vdeh3'};
model=addReaction(model,'az','Ax -> Zx'); cat(length(model.rxns))={'vdeh3'};
model=addReaction(model,'za','Zx -> Ax');
model=addReaction(model,'av','Ax -> Vx');
model=addReaction(model,'zxpsbs','Zx + psbsh6 <=> psbsh6zx');
model=addReaction(model,'q','P680e -> P680'); cat(length(model.rxns))={'psbsh6zx'};


% State transition
model=addReaction(model,'stn7_ac','stn7 + 2 PQH2 <=> stn7a + 2 PQ');
model=addReaction(model,'pph1_ac','pph1 + 2 PQ <=> pph1a + 2 PQH2');
model=addReaction(model,'chl12','chl1_m -> chl2_m'); cat(length(model.rxns))={'pph1a'};
model=addReaction(model,'chl21','chl2_m -> chl1_m'); cat(length(model.rxns))={'stn7a'};


% Exchange, sink and source

model=addReaction(model,'kea3','H[l] + k[l] -> k[s]'); 
model=addReaction(model,'o2ex','O2[l] <=> O2[s]'); 
model=addReaction(model,'o2ex2','O2[s] <=> O2[c]'); 
model=addReaction(model,'o2ext',' O2[c] <=> '); 

model=addReaction(model,'co2ex','CO2[c] <=> CO2[s]'); 
model=addReaction(model,'co2ex2',' <=> CO2[c]'); 

model=addReaction(model,'ho2sl','H2O[l] <=> H2O[s]'); 
model=addReaction(model,'ho2sc','H2O[c] <=> H2O[s]'); 
model=addReaction(model,'h2osink','H2O[c] <=> '); 

model=addReaction(model,'mgex','mg2p[s] <=> '); 

model=addReaction(model,'Eex','E <=> '); 

model=addReaction(model,'Hsc','H[c] <=> H[s]'); 
model=addReaction(model,'Pisc','Pi[c] <=> Pi[s]');
model=addReaction(model,'ATPexc',' <=> ATP[s]');

model=addReaction(model,'sink_SUC','SUC <=> ');


% for last reaction: extend the vector cat to the number of reactions
cat(length(model.rxns))={''}; 

model= rmfield(model,'b');
model.rev=zeros(length(model.rxns),1);
model.rev(model.lb<0)=1;