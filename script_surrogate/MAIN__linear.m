
load('initial_data_new.mat')

%% PARAMETERIZATION

% STEP3: surrogate model

% Parsing the structure data
model_d.leftout=5;
interp=0;
tn=[];

model_d.nm=length(model_d.mets); % number of metabolites
model_d.nr=length(model_d.rxns); % number of reactions
model_d.ncond=5-1; % number of conditions
model_d.conditions=setdiff(1:5,model_d.leftout);
model_d.k0=model_d.ncond*model_d.nm; % position where k starts

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
    if cond~=model_d.leftout
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


%% Constraints for surrogate optimization

% 1) linear inequality constraints in form of A*x <= b

% Inequality constraints
% 1) sum of enzyme and enzyme complexes to 1
% e_j,1 + e_j,2 + ... + e_j,n =e_j_total<=1;


A1=zeros(size(enzymes.name,1),model_d.nm);
for enz=1:size(enzymes.name,1)
    A1(enz,model_d.cxe(:,enz)~=0)=1; 
end

Ar = repmat(A1, 1, 4);                                  
Ac = mat2cell(Ar, size(A1,1), repmat(size(A1,2),1,4));
A1 = sparse(blkdiag(Ac{:}));

b1=ones(size(A1,1),1);


A=A1;
b=b1;
% 2) linear equality constraints in form of A*x = b

Aeq=[];
beq=[];


%% Lower and upper bounds
lb=[repmat(model_d.metslb,4,1);klbub(1)]; 
ub=[repmat(model_d.metsub,4,1);klbub(2)];


%% Objective function including non-linear constraints
objconstr = @(x)SURR_obj_const(x,model_d);


% non linear constraint:
% Applied for the initial time points (steady state)
% 1) Steady state
% 2) circle property of reaction rate constants k at steady state


init_sol=[reshape(init_meta(:,model_d.conditions),[],1);init_k];
options = optimoptions('surrogateopt','MaxFunctionEvaluations',20,'PlotFcn','surrogateoptplot','UseParallel',true,'InitialPoints',init_sol,'MinSampleDistance',1);

options.PlotFcn = 'surrogateoptplot';
options.UseParallel = true;

intcon=[];
problem = struct('objective',objconstr,'lb',lb,'ub',ub, 'Aineq', A, 'bineq',b,...
    'options',options,'solver','surrogateopt');

[x,fval,exitflag,output,trials] = surrogateopt(problem);
