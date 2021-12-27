function z = obj_function(x,model_d)

%% OBJECTIVE FUNCTION
%nt=model_d.nt;
nm=model_d.nm;

%timept=model_d.timept;
k0=model_d.k0;


log_c=x(1:k0);
k=10.^(x(k0+1:end));

%%%Integration of metabolite concentration at next timepoint
concentrations=cell(model_d.ncond,1);
%zeros(nm,nt,model_d.ncond); % 3D matrix of metabolites: rows(metabolites), columns(time points), width(experiments)
% ! The initial values are in log scale for metabolites and ks
% So it need to be transformed back later
log_initial_conc=reshape(log_c,nm,model_d.ncond);

% dimension: mets x time x cond

% ! The concentration matrices are in the original scale
linear_x=[];

for cond=1:model_d.ncond
    concentrations{cond}=zeros(nm,size(model_d.c3d_meas{cond},2)); % mets x time points
    
    backup=log_initial_conc(model_d.enz_ind,cond);
    initial_conc=10.^(log_initial_conc(:,cond));
    initial_conc(model_d.enz_ind,:)=backup;
    
    concentrations{cond}(:,1)=initial_conc;
    linear_x=[linear_x,initial_conc];
end   


% the time interval for the interpolation was in hours, then: 
% htomin=60;  because the ks are estimated in minutes
htomin=1; % because the ks are estimated in hours

%%
% for every condition

for cond=1:model_d.ncond
    sol=[];
    tspan=model_d.irra_ode{cond}(:,2)*htomin;
    x0=concentrations{cond}(model_d.sim_ind,1);
    
    startT=datetime('now');
    overfunc=@(sim_t,sim_x)stopfunc(sim_t,sim_x,startT);
    options=odeset('NonNegative',1:length(model_d.sim_ind),'RelTol',1e-3,'AbsTol',1e-3,'Events',overfunc,'InitialStep',0.001);

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

z=0;
concentrations_pred=cell(model_d.ncond,1);
for cond=1:model_d.ncond
    concentrations_pred{cond}=concentrations{cond}(model_d.model_ind,:);

    % Calculate the sum of square difference
    sqsum=(reshape(model_d.c3d_meas{cond},[],1)-reshape(concentrations_pred{cond},[],1)).^2./maxsd;
    z=z+nansum(sqsum);
end

