function [res]=get_DAE_new(t,sim_x,sim_dx,kinetic_param,model_d,cond,dae_ind)

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

dxdt_full(3)=x_full(2)+x_full(3)-1;
dxdt_full(5)=x_full(4)+x_full(5)-1;
dxdt_full(7)=x_full(6)+x_full(7)-1;
dxdt_full(9)=x_full(8)+x_full(9)-1;
dxdt_full(11)=x_full(10)+x_full(11)-1;
dxdt_full(13)=x_full(12)+x_full(13)-1;
% dxdt_full(18)=sum(x_full([14,15,18]))-1; % P680
dxdt_full(39)=sum(x_full([27:35,39]))-1; % S0
dxdt_full(46)=sum(x_full([40,41,43,46]))-1; % ISP
dxdt_full(52)=sum(x_full(49:52))-1; % cytochrome
%  dxdt_full(59)=sum(x_full([55,56,59]))-1; % P700

for enz=1:4% size(model_d.enzymes_name,1)
    ind=find(model_d.cxe(:,enz)~=0);
    dxdt_full(ind(end))=sum(x_full(ind))-1; 
end

% only the simulated metabolite are returned as output
dxdt=dxdt_full(sim_ind);


sim_dx(dae_ind)=0;
res=sim_dx-dxdt;



