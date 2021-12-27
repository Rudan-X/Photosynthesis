function [init_flux]=get_init_flux(model_d,light_5)


init_flux=zeros(length(model_d.rxns),5);
%% Geometric FBA

for cond=1:5
    model_d.c=zeros(length(model_d.rxns),1);
    
    model_d.ub=1e5*ones(size(model_d.S,2),1);
    %model_d.lb=zeros(size(model_d.S,2),1);
    model_d.lb=1e-3*ones(size(model_d.S,2),1);
    
    model_d.lb(1)=0.999*light_5(cond);
    model_d.ub(1)=1.001*light_5(cond);

    [~,ind]=ismember('carb1',model_d.rxns);
    %[~,ind]=ismember('carb1',model_d.rxns);
    model_d.c(ind)=1;
    
    [~,ind]=ismember('ss1_f',model_d.rxns);
    % [~,ind]=ismember('ss',model_d.rxns);
    model_d.c(ind)=1;
    
    % [~,ind]=ismember('ss_b',model_d.rxns);
    % [~,ind]=ismember('ss1_b',model_d.rxns);
    % model_d.c(ind)=-1;

%     [~,ind]=ismember('q1_f',model_d.rxns);
%     model_d.c(ind)=-1;
        
    init_flux(:,cond) = geometricFBA(model_d);
end
