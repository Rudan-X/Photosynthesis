function [dxdt]=get_ode(t,sim_x,kinetic_param,model_d,cond)

fix_ind=model_d.fix_ind;
sim_ind=model_d.sim_ind;
ind=intersect(find(model_d.irra_ode{cond}(:,1)==0),find(model_d.irra_ode{cond}(:,2)>0));
if isempty(ind)
    fix_x=model_d.constant(fix_ind,1);
else
    ind=ind(1);
    night_t=model_d.irra_ode{cond}(ind,2);
    if t<night_t
        fix_x=model_d.constant(fix_ind,1);
    else % when there is no light
        fix_x=model_d.constant(fix_ind,2);
    end
end


% reconstruct the full species matrix:
x_full=zeros(length(sim_x)+length(fix_ind),1);
x_full(sim_ind)=sim_x;
x_full(setdiff(1:length(x_full),sim_ind))=fix_x;


v=kinetic_param.*prod(repmat(x_full,1,length(kinetic_param)).^(abs(model_d.S.*(model_d.S<0))))';
ind=find(model_d.irra_ode{cond}(:,2)<=t,1,'last');
v(1)=model_d.irra_ode{cond}(ind,1);

dxdt_full=model_d.S*v;
% only the simulated metabolite are returned as output
dxdt=dxdt_full(sim_ind);

