


tspan=model_d.irra_ode{cond}(:,2)*htomin;
x0=concentrations{cond}(model_d.sim_ind,1);
dx0=zeros(length(x0),1);

startT=datetime('now');
overfunc=@(sim_t,sim_x)stopfunc(sim_t,sim_x,startT);

toadd=[];
for enz=1:4 %size(model_d.enzymes_name,1)
    ind=find(model_d.cxe(:,enz)~=0);
    toadd=[toadd,ind(end)];
end


ind=[3,5,7,9,11,13,39,46,52,toadd];
[~,DAE_ind]=ismember(ind,model_d.sim_ind);
jacob=diag(ones(length(model_d.sim_ind),1));
jacob(DAE_ind,:)=0;
options=odeset('RelTol',1e-1,'AbsTol',1e-1,'Jacobian',{[],jacob});
sol = ode15i(@(sim_t,sim_x,sim_dx)get_DAE_new(sim_t,sim_x,sim_dx,k,model_d,cond,DAE_ind),tspan,x0,dx0,options);



%%%%%%%%%%%% ODE15S %%%%%%%%%%%%%%
M = diag(ones(length(model_d.sim_ind),1));
ind=[3,5,7,9,11,13]; % index from model_d.mets
[~,ind]=ismember(ind,model_d.sim_ind);
M(ind,:)=0;
options=odeset('Mass',M,'RelTol',1e-2,'AbsTol',1e-2,'Events',overfunc,'InitialStep',0.001);
sol=ode15s(@(sim_t,sim_x)get_ode_new(sim_t,sim_x,k,model_d,cond),tspan,x0,options); 



options=odeset('RelTol',1e-2,'AbsTol',1e-2,'Events',overfunc,'InitialStep',0.001);
sol=ode15s(@(sim_t,sim_x)get_ode(sim_t,sim_x,k,model_d,cond),tspan,x0,options); 


sol=ode15s(@(sim_t,sim_x)get_ode_old(sim_t,sim_x,k,model_d,cond),tspan,x0,options); 

