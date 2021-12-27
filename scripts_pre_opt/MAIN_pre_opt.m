% changeCobraSolver('gurobi');

%% CREATING THE MODEL
Enzymes() % list of enzymes
% original reactions
build_model()

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
irra_ini=zeros(length(model_d.irradiance),2);
for cond=1:length(model_d.irradiance)
    ind=find(model_d.irradiance{cond}(:,1)~=0,1);
    irra_ini(cond,1)=model_d.irradiance{cond}(ind,1); % flux value
    irra_ini(cond,2)=ind; % index of the time point in the matrix
    irra_ini(cond,3)=model_d.irradiance{cond}(ind,2); % time point in hour
end
%irra_ini=[irradiation, index in the vector, time point], all at the
%starting time point


[init_flux]=get_init_flux(model_d,irra_ini(:,1));
% returned flux units: mili mol / L / h: mM/h
%% Lower and upper bounds
% general bounds for metabolites
metslb=1e-8;
metsub=1e8;
model_d.metslb=metslb*ones(length(model_d.mets),1);
model_d.metsub=metsub*ones(length(model_d.mets),1);

% the upper bounds for enzymes and enzyme complexes, and
% also their non-catalytic forms, and those whose reactions are manually 
% decomposed are set to 1

is_enz_comp=sum(model_d.cxe,2); 
model_d.metsub(is_enz_comp==1)=1;


% lower und upper bounds for k and epsilon
klbub=[1e-6, 1e6]; 
epslbub=[1e-10,1e10];
[init_meta,init_k,eps]=get_initials(model_d,init_flux,irra_ini(:,3),klbub,epslbub);

save('init_data.mat')