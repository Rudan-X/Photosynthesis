load('initial_data_new.mat')

metslb=1e-8;
metsub=1e8;
model_d.metslb=metslb*ones(length(model_d.mets),1);
model_d.metsub=metsub*ones(length(model_d.mets),1);

enz_ind=find(sum(model_d.cxe,2)==1); % enzymes and enzyme complexes, enz_ind are all found in sim_ind

model_d.metsub(enz_ind)=1;

%%  Parsing the structure data
model_d.leftout=5;
interp=0;
tn=10;

model_d.nm=length(model_d.mets); % number of metabolites
model_d.nr=length(model_d.rxns); % number of reactions
model_d.ncond=5-1; % number of conditions
model_d.conditions=setdiff(1:5,model_d.leftout);
model_d.k0=model_d.ncond*model_d.nm; % position where k starts
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
        if isempty(tn)
            index=irra_ini(cond,2):size(data_cell{cond,1},1);
        else
            index=irra_ini(cond,2):(irra_ini(cond,2)+tn);
        end
        model_d.c3d_meas{i}=data_cell{cond,1}(index,table_ind)';
        model_d.sd3d{i}=data_cell{cond,2}(index,table_ind)';
        model_d.irra_ode{cond}=model_d.irradiance{cond}(index,:);
    end
end
