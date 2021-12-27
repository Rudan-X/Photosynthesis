function [meta,rates,epsi]=get_initials(model_d,fluxes,time,klbub,epslbub)
%% Variable vector: [meta_concentration (condition 1:num_cond), k(unique), epsilon (condition 1:num_cond)]
num_cond = 5;
[m,n]=size(model_d.S);
nvar=num_cond*m+n+num_cond*n; % number of variables

%% First equality constraint
% v=n * c^k is reformulated into:
% log(v)=log(k)+n*log(c)+ep;

Aeq1=zeros(num_cond*n,nvar);
beq1=zeros(num_cond*n,1);

for cond=1:num_cond
    k0=num_cond*m + 1;
    kend=num_cond*m + n;
    row=(cond-1)*n+1:cond*n;
    col=(cond-1)*m+1:cond*m;
    eprange=num_cond*m+n+(cond-1)*n+1:num_cond*m+n+(cond-1)*n+n;

    temp=model_d.S;
    %get only the positive entries
    temp(temp>0)=0;
    temp=abs(temp);

    Aeq1(row,col)=temp';
    Aeq1(row,k0:kend)=eye(n);
    Aeq1(row,eprange)=eye(n);

    temp=log(fluxes(:,cond));
    beq1(row)=temp;

end

%% Second equality constraints
% log(k1*k2*k3/k-1/k-2/k-3)=log(k1)+log(k2)+log(k3)-log(k-1)-log(k-2)-log(k-3)

Aeq2=zeros(max(model_d.cata),nvar);
beq2=log(ones(max(model_d.cata),1));

for i=1:max(model_d.cata)
    range=num_cond*m +1:num_cond*m + n;
    if sum(model_d.cata==i)==sum(model_d.cata==-i) %only for those with reversible cycles
        %check
        %display(model_d.rxns(find(model_d.cata==i)))
        Aeq2(i,range(model_d.cata==i))=1;
        Aeq2(i,range(model_d.cata==-i))=-1;
    end
end

%%
prob=struct();

%% General bounds
% metabolite concentration (Already stored in model_d.metsub and .metslb): 

prob.lb=[log(repmat(model_d.metslb,num_cond,1));log(klbub(1))*ones(n,1);log(epslbub(1))*ones(num_cond*n,1)]; 
prob.ub=[log(repmat(model_d.metsub,num_cond,1));log(klbub(2))*ones(n,1);log(epslbub(2))*ones(num_cond*n,1)];

%% 3. experimental metabolite => upper and lower bounds for each of the num_cond conditions
full_meta=model_d.exp_meta;

[~,ind]=ismember(full_meta{6,1},string(model_d.mets));
table_ind=find(ind~=0);
mets_ind=ind(table_ind); % indices in the model

for cond=1:num_cond 
    time_ind=find(full_meta{cond,1}(:,end)==time(cond));
    time_ind=time_ind(1);
    var_ind=(cond-1)*m+mets_ind;
    temp=(1-1e-7)*full_meta{cond,1}(time_ind,table_ind);

    % check if some experimental data are missing, if it is missing, we do
    % not consider it for bounds (reset to default value)
    temp(isnan(temp))=0;
    prob.lb(var_ind)=log(temp);
    temp=(1+1e-7)*full_meta{cond,1}(time_ind,table_ind);
    
    temp(isnan(temp))=max(model_d.metsub);
    temp(temp==0)=max(model_d.metsub);
    prob.ub(var_ind)=log(temp);
    
end



%% Objective min sum((epsilon)^2)
prob.obj=zeros(nvar,1);
tempQ=zeros(nvar);
tempQ(num_cond*m+n+1:end,num_cond*m+n+1:end)=eye(num_cond*n);


prob.Q=sparse(tempQ);

Aeq=[Aeq1;Aeq2];
beq=[beq1;beq2];

Aineq=[]; bineq=[];
prob.A=[sparse(Aeq); sparse(Aineq)];
prob.rhs = full([beq;bineq]);
prob.vtype = repmat('C', nvar, 1);
prob.sense = [repmat('=',size(Aeq,1),1);repmat('<',size(Aineq,1),1)];


sol=gurobi(prob);

meta=exp(reshape(sol.x(1:num_cond*m),m,num_cond));
rates=exp(reshape(sol.x(num_cond*m+1:num_cond*m+n),n,1));
epsi=exp(reshape(sol.x(num_cond*m+n+1:end),n,num_cond));

% %% Variability of k
% prob2=prob;
% prob2=rmfield(prob2,'Q');
% prob2.quadcon.Qc = sparse(tempQ);
% prob2.quadcon.q  = zeros(nvar,1); 
% prob2.quadcon.rhs = sol.objval;
% prob2.quadcon.sense = repmat('<',length(prob2.quadcon.rhs),1);
% prob2.quadcon.name = 'rot_cone';
% 
% minmax_k=zeros(n,2);
% 
% for reac=21:n
%     objf=zeros(nvar,1);
%     objf(5*m+reac)=1;
%     prob2.obj=objf;
%     sol_min=gurobi(prob2);
%     minmax_k(reac,1)=exp(sol_min.objval);
%     objf(5*m+reac)=-1;
%     prob2.obj=objf;
%     sol_max=gurobi(prob2);    
%     minmax_k(reac,2)=exp(-sol_max.objval);
% end
end