function [Fval,g]=MEI_obj_const_log(x,model_d)

nm=model_d.nm;

%timept=model_d.timept;
k0=model_d.k0;

k=10.^(x(k0+1:end));
log_initial_conc=x(1:k0);
%%%Integration of metabolite concentration at next timepoint

%%%%%%%%%%%%%%%%%%%%%%%%
concentrations=cell(model_d.ncond,1);
%zeros(nm,nt,model_d.ncond); % 3D matrix of metabolites: rows(metabolites), columns(time points), width(experiments)
% ! The initial values are in log scale for metabolites and ks
% So it need to be transformed back later


% dimension: mets x time x cond

% ! The concentration matrices are in the original scale
linear_x=[];
for cond=1:model_d.ncond
    concentrations{cond}=zeros(nm,size(model_d.c3d_meas{cond},2)); % mets x time points
    
    backup=log_initial_conc(model_d.enz_ind);
    initial_conc=10.^(log_initial_conc);
    initial_conc(model_d.enz_ind)=backup;
    
    concentrations{cond}(:,1)=initial_conc;
    linear_x=[linear_x,initial_conc];
end

%%%%%%%%%%%%%%%%%%%%%%%%


% the time interval for the interpolation was in hours, then: 
% htomin=60;  because the ks are estimated in minutes
htomin=1; % because the ks are estimated in hours

%%
% for every condition
CBC_met={'PGA[s]','ATP[s]','ADP[s]','DpGA[s]', 'FBP[s]', 'H2O[s]',  'F6P[s]','Pi[s]',...
    'F6P[s]','NADPH[s]', 'NADPp[s]','GAP[s]', 'DHAP[s]', 'E4P','X5P', 'SBP', 'S7P',...
    'R5P', 'X5P', 'RU5P', 'RuBP',  'G6P[s]', 'G1P[s]','ADPG', 'glucan', 'STA'};


DAE_ind=contains(string(model_d.mets),string(CBC_met));
[~,DAE_ind]=ismember(DAE_ind,model_d.sim_ind);

DAE_ind=setdiff(1:length(model_d.sim_ind),DAE_ind);
DAE_ind=DAE_ind(6);
jacob=diag(ones(length(model_d.sim_ind),1));
jacob(DAE_ind,:)=0;

for cond=1:4
    
    fprintf(' Condition: %d\n',cond);
    
    tspan=model_d.irra_ode{cond}(:,2)*htomin;
    x0=concentrations{cond}(model_d.sim_ind,1);
    dx0=zeros(length(x0),1);
    
    startT=datetime('now');
    overfunc=@(sim_t,sim_x)stopfunc(sim_t,sim_x,startT);
 
    options=odeset('RelTol',1e-1,'AbsTol',1e-1,'Jacobian',{[],jacob});
    sol = ode15i(@(sim_t,sim_x,sim_dx)get_DAE(sim_t,sim_x,sim_dx,k,model_d,cond,DAE_ind),tspan,x0,dx0,options);

    options=odeset('RelTol',1e-2,'AbsTol',1e-2,'Events',overfunc,'InitialStep',0.001);
    sol=ode15s(@(sim_t,sim_x)get_ode(sim_t,sim_x,k,model_d,cond),tspan,x0,options); 

    % if sol.ie is empty, the integration is completed correctly. If
    % sol.ie==1, the integration is terminated without completing the
    % whole integration time interval.
    if isempty(sol)
        concentrations{cond}(model_d.sim_ind,ind+1:end)=0;
    else
        maxtime=max(sol.x);
        ind=find(tspan<=maxtime,1,'last');
        concentrations{cond}(model_d.sim_ind,1:ind)=deval(sol,tspan(1:ind));
        concentrations{cond}(model_d.sim_ind,ind+1:end)=0;
    end
    
end        
%%
% check the metabolites with measured concentrations

maxsd=0;
for cond=1:model_d.ncond
    maxsd=max([maxsd,max(model_d.sd3d{cond},[],"all")]); % maximum standard deviation
end

Fval=0;

for cond=1:model_d.ncond
    concentrations_pred=concentrations{cond}(model_d.model_ind,2:end);

    % Calculate the sum of square difference
    sqsum=(reshape(model_d.c3d_meas{cond}(:,2:end),[],1)-reshape(concentrations_pred,[],1)).^2./maxsd;
    Fval=Fval+nansum(sqsum);
end

% R=(yout-yexp);
% R=reshape(R,numel(R),1);



%%%%%%%%%%%%%%% EQUALITY CONSTRAINT
%% Constraints for optimization

%%%%%%%%%% LINEAR CONSTRAINTS%%%%%%%%%%%%%%%
% linear inequality constraints in form of A*x <= b
% 1) sum of enzyme and enzyme complexes to 1
% e_j,1 + e_j,2 + ... + e_j,n =e_j_total<=1;

A1=zeros(size(model_d.enzymes_name,1),model_d.nm);
for enz=1:size(model_d.enzymes_name,1)
    A1(enz,model_d.cxe(:,enz)~=0)=1; 
end

% Ar = repmat(A1, 1, model_d.ncond);                                  
% Ac = mat2cell(Ar, size(A1,1), repmat(size(A1,2),1,model_d.ncond));
% A1 = sparse(blkdiag(Ac{:}));

b1=ones(size(A1,1),1);

A=[A1,zeros(size(A1,1),model_d.nr)];

g1=A*x'-b1;

%%%%%%%%%% NONLINEAR CONSTRAINTS%%%%%%%%%%%%%%%

% 1) steady state applied only for the INITIAL TIME POINT

% Ar = repmat(model_d.S, 1, model_d.ncond);                                   % Repeat Matrix
% Ac = mat2cell(Ar, size(model_d.S,1), repmat(size(model_d.S,2),1,model_d.ncond));    % Create Cell Array Of Orignal Repeated Matrix
% S = sparse(blkdiag(Ac{:}));

S=model_d.S;
subst=S.* (S<0);

k0=model_d.k0;
k=10.^(x(k0+1:end))';

meta=linear_x(:,1);
k2=k; 
cmap=repmat(reshape(meta,[],1),1,length(k2));
flux=k2.* prod(cmap.^(abs(subst)))';

geq1=S*flux; 

% 2) kinetic parameters relationship inside circles at steady state
geq2=zeros(max(model_d.cata),1);
for i=1:max(model_d.cata)
    if sum(model_d.cata==i)==sum(model_d.cata==-i) %only for those with reversible cycles
        geq2(i)=prod(k(model_d.cata==i))/prod(k(model_d.cata==-i))-1;
    end
end

g=[g1;geq1;geq2];

return

%***************************************************
%Function of the dynamic system



   
     
   
 
 
 