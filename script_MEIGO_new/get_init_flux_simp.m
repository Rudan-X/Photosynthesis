function [init_flux]=get_init_flux_simp(model_d,light)


%% Geometric FBA

model_d.c=zeros(length(model_d.rxns),1);

model_d.ub=1e5*ones(size(model_d.S,2),1);
%model_d.lb=zeros(size(model_d.S,2),1);
model_d.lb=1e-3*ones(size(model_d.S,2),1);

model_d.lb(1)=0.999*light;
model_d.ub(1)=1.001*light;

[~,ind]=ismember('carb1',model_d.rxns);
% [~,ind]=ismember('carb1',model_d.rxns);
model_d.c(ind)=1;

% [~,ind]=ismember('ss1_f',model_d.rxns);
[~,ind]=ismember('ss',model_d.rxns);
model_d.c(ind)=1;


% [~,ind]=ismember('ss1_b',model_d.rxns);
% model_d.c(ind)=-1;

%     [~,ind]=ismember('q1_f',model_d.rxns);
%     model_d.c(ind)=-1;

init_flux = geometricFBA(model_d);