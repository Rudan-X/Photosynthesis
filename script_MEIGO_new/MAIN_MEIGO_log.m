pre_optimization()

cd ..
cd meigo64-2020/MEIGO
install_MEIGO
cd ..
cd ..
cd script_MEIGO/

addpath 'work/xu2/meigo64-2020\MEIGO\eSS'
% addpath 'C:\Users\Ruru\Documents\MATLAB\meigo64-2020\MEIGO\eSS'
%% Lower and upper bounds
% log scale:

model_d.metslb(model_d.model_ind)=(1-0.001)*init_meta(model_d.model_ind);
backup=model_d.metslb(enz_ind);
log_lb=log10(model_d.metslb);
log_lb(enz_ind)=backup;

model_d.metsub(model_d.model_ind)=(1+0.001)*init_meta(model_d.model_ind);
backup=model_d.metsub(enz_ind);
log_ub=log10(model_d.metsub);
log_ub(enz_ind)=backup;

lb=[log_lb;repmat(log10(klbub(1)),model_d.nr,1)]; 
ub=[log_ub;repmat(log10(klbub(2)),model_d.nr,1)];


% Objective function including non-linear constraints

backup=init_meta(enz_ind);
log_init_x=log10(init_meta);
log_init_x(enz_ind)=backup;
init_sol_log=[log_init_x;log10(init_k)];

% init_sol_linear=[reshape(init_meta(:,model_d.conditions),[],1);init_k];
% 
% check=load('patternsearch_result_3days.mat');
% init_sol_log=check.x;

%%
problem.f='MEI_obj_const_log';           
problem.x_L=lb;
problem.x_U=ub;

problem.x_0=init_sol_log;
opts.maxeval=20; 
% opts.local.solver='nl2sol';
opts.local.solver='fmincon';
% opts.local.solver='lsqnonlin';
opts.inter_save=1;


Results=MEIGO(problem,opts,'ESS',model_d);
save('res_MEIGO.mat','Results')
