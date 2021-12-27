% load('initial_data_new.mat')

changeCobraSolver('gurobi');

%% CREATING THE MODEL
Enzymes() % list of enzymes
% Enzymes_simp()
% original reactions
build_model()

% build_model_simp()

% decompose the model
[model_d,model0]=decomp_model(model,enzymes,cat);

% add ratio between carboxylation and photorespiration
model_d = addRatioReaction (model_d, {'carb1', 'resp1'}, [1; 1.5]);

% add ratio between starch and sucrose synthesis
model_d = addRatioReaction (model_d, {'ss1', 'spase1'}, [1; 1]);

% model_d = addRatioReaction (model_d, {'ss', 'spase'}, [1; 1]);

% convert the reversible model into irreversible one
model_d=convertToIrreversible(model_d);
model_d.cata(find(contains(string(model_d.rxns),'_b'),1):end)=-model_d.cata(find(contains(string(model_d.rxns),'_b'),1):end);

% chech if there is any blocked reaction
block=findBlockedReaction(model_d)



%%
% STEP1: obtain initial steady state flux distribution
% Some fluxes are constrained with amount/time: inside the function
% get_init_flux, the units of fluxes are changed to mol/L/s
% In case other units are wanted:
% number of seconds per the time unit requested:
% number of mol per the amount unit requested:

second_per_time_unit=3600; % 3600 second per hour: returned time units=hour
mol_per_amount_unit=1e-3; % 0.001 mili mol per mol: returned amount units= mili mol



% Load experimental data
[model_d]=load_experiment_data(model_d,mol_per_amount_unit); % store the info in model_d.exp_meta
% model_d.exp_meta=[concentration, std dev]

%Record the metabolites with fixed concentrations and reactions with fixed fluxes
[model_d,metacte]=load_ODE_cte(model_d,mol_per_amount_unit,second_per_time_unit); % model_d.constant=[day,night]


% find the first time point the irradiance is non-zero

% irra_ini=zeros(length(model_d.irradiance),2);
% for cond=1:length(model_d.irradiance)
%     ind=find(model_d.irradiance{cond}(:,1)~=0,1);
%     irra_ini(cond,1)=model_d.irradiance{cond}(ind,1); % flux value
%     irra_ini(cond,2)=ind; % index of the time point in the matrix
%     irra_ini(cond,3)=model_d.irradiance{cond}(ind,2); % time point in hour
% end

%irra_ini=[irradiation, index in the vector, time point], all at the
%starting time point

mean=zeros(5,1);
for cond=1:length(model_d.irradiance)
    ind=find(model_d.irradiance{cond}(:,2)==0,1);
    mean(cond)=model_d.irradiance{cond}(ind,1);
end

cond=1;
ind=find(model_d.irradiance{cond}(:,2)==0,1);

[init_flux]=get_init_flux(model_d,model_d.irradiance{cond}(ind,1)); %irra_ini(1,1)


% [init_flux]=get_init_flux_simp(model_d,model_d.irradiance{cond}(ind,1)); %irra_ini(1,1)

% returned flux units: mili mol / L / h: mM/h
%% Lower and upper bounds
% general bounds for metabolites
metslb=1e-6;
metsub=1e6;
model_d.metslb=metslb*ones(length(model_d.mets),1);
model_d.metsub=metsub*ones(length(model_d.mets),1);

% the upper bounds for enzymes and enzyme complexes, and
% also their non-catalytic forms, and those whose reactions are manually 
% decomposed are set to 1

enz_ind=find(sum(model_d.cxe,2)==1); % enzymes and enzyme complexes, enz_ind are all found in sim_ind

model_d.metsub(enz_ind)=1;


% lower und upper bounds for k and epsilon
klbub=[1e-6, 1e6]; 
epslbub=[1e-10,1e10];

mat=zeros(5,22);
for cond=1:5
    mat(cond,:)=model_d.exp_meta{cond,1}(1,1:22);
end
mean_meta0=nanmean(mat,1);
[init_meta,init_k,eps]=get_initials(model_d,init_flux,mean_meta0,klbub,epslbub);


v=init_k.*prod(repmat(init_meta,1,length(init_k)).^(abs(model_d.S.*(model_d.S<0))))';
% comp=[init_flux,v,eps,string(model_d.rxns)];

%%  Parsing the structure data
model_d.leftout=5;
interp=0;
tf=16;

model_d.nm=length(model_d.mets); % number of metabolites
model_d.nr=length(model_d.rxns); % number of reactions
model_d.ncond=5-1; % number of conditions
model_d.conditions=setdiff(1:5,model_d.leftout);
model_d.k0=model_d.nm; % position where k starts
model_d.enz_ind=enz_ind;
model_d.enzymes_name=enzymes.name;

model_d.fix_ind=find(model_d.constant(:,1)~=0); % fixed metabolites
model_d.sim_ind=find(model_d.constant(:,1)==0); % simulated metabolites

[~,enz_in_sim]=ismember(enz_ind,model_d.sim_ind);
model_d.enz_in_sim=enz_in_sim;
% number of variables
nvars=model_d.k0+model_d.nr;


% Time points evaluated:
% interpolation option included
if interp==1
    data_cell=model_d.exp_meta_interp;
elseif interp==0
    data_cell=model_d.exp_meta;
end



%%
[~,ind]=ismember(data_cell{6,1},string(model_d.mets));
table_ind=find(ind~=0);
model_d.model_ind=ind(table_ind);

model_d.c3d_meas=cell(model_d.ncond,1);
model_d.sd3d=cell(model_d.ncond,1);


i=0;
model_d.irra_ode=cell(model_d.ncond,1);
for cond=1:model_d.ncond+1
    % [~,tm]=ismember(timept,data_cell{cond,1}(:,end));
    if sum(cond==model_d.conditions)~=0
        i=i+1;
        if isempty(tf)
            index=1:size(data_cell{cond,1},1);
        else
            k=find(data_cell{cond,1}(:,end)==tf);
            index=1:k;
        end
        model_d.c3d_meas{i}=data_cell{cond,1}(index,table_ind)';
        model_d.sd3d{i}=data_cell{cond,2}(index,table_ind)';
        model_d.irra_ode{cond}=model_d.irradiance{cond}(index,:);
    end
end
