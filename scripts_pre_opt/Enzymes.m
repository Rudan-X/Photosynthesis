enzymes=struct();
enzymes.name={'atpsac';...
    'fnr';...
    'ftr';...
    'fbpase';...
    'fbpasen';...
    'sbpase';...
    'sbpasen';...
    'prk';...
    'prkn';...
    'gapdh_s';...
    'pgak_s';...
    'tpi_s';...
    'fbpaldo_s';...
    'transk';...
    'rpi';...
    'rpe';...
    'gpi_s';...
    'pgm_s';...
    'adpgpy';...
    'starsyn';...
    'agb';...
    'spase';
    'f26bpase';...
    'pfk2';...
    'ndpk';
    'pgak_c';...
    'gapdh_c';...
    'tpi_c';...
    'fbpaldo_c';...
    'fbpase_c';...
    'gpi_c';...
    'pgm_c';...
    'ugpase';...
    'sps';...
    'vdeh3';...
    'psbsh6zx';...
    'stn7a';...
    'pph1a'};

enzymes.EC={'3.6.3.14';...
    '1.18.1.2';...
    '1.8.7.2';...
    '3.1.3.11';...
    '3.1.3.11';...
    '3.1.3.37';...
    '3.1.3.37';...
    '2.7.1.19';...
    '2.7.1.19';...
    '1.2.1.13';... % ping pong
    '2.7.2.3';...
    '5.3.1.1';...
    '4.1.2.13';... % Keto acid
    '2.2.1.1';...
    '5.3.1.6';...
    '5.1.3.1';...
    '5.3.1.9';...
    '5.4.2.2';...
    '2.7.7.27';... % Phosphorous-containing groups
    '2.4.1.21';... %  Glycosyl groups
    '2.4.1.1';... %  Glycosyl groups
    '3.1.3.24';
    '3.1.3.46';...
    '2.7.1.105';...
    '2.7.4.6';
    '2.7.2.3';...
    '1.2.1.13';... %  Aldehyde or oxo groups
    '5.3.1.1';...
    '4.1.2.13';...
    '3.1.3.11';... % Ester bonds
    '5.3.1.9';...
    '5.4.2.2';...
    '2.7.7.9';...
    '2.4.1.14';...
    '1.14.15.21';...
    '';...
    '2.7.11.1';...
    '3.1.3.16'};
enzymes.noncata={'P680';'P680e';'pheo';'pheon';'P680p';'QA';'QAn';'QB';'QBn';'QBn2';'s0';'s0p';'s1';'s1p';'s2';'s2p';'s3';'s3p';'s4';'s4p';...
    'ISPox';'ISPoxPQH2';'PQsemi';'ISPHred';'cytf';'cytfn';'ISPHp';'cytbl';'cytbln';'cytbh';'cytbhn';'P700';'P700e';'A0';'A0n';'P700p';...
    'atps';'Rcai';'Rcaa';'ER';'ERRcaa';'ERcaa';'E';'ECR';'ECRRcaa';'ECRcaa';'EC';'ECM';'ECMR';'ECMRC';'RPGA2H';'RPGA2';'RPGA';...
    'ECMRO2';'PGAPGLRubiH2';'PGAPGLRubi';'PGAPGL'; 'psbs';'psbsh6';'VDE';'stn7';'pph1'};

% 1: oxidoreductase
% 2: Transferase
% 3: Hydrolase: hydrolytic cleavage
% 4: Lyases: remove bonds(aldolase, decarboxylase)
% 5: Isomerase: geometric or structural changes in 1 molecule
% 6: synthase: catalyze the junction of 2 molecules
% 7: translocase: movement of ions or molecules across membrane
