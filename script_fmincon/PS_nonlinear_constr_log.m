function [c,ceq] = PS_nonlinear_constr_log(x,model_d)

c = [];

% 1) steady state applied only for the INITIAL TIME POINT

Ar = repmat(model_d.S, 1, 4);                                   % Repeat Matrix
Ac = mat2cell(Ar, size(model_d.S,1), repmat(size(model_d.S,2),1,4));    % Create Cell Array Of Orignal Repeated Matrix
S = sparse(blkdiag(Ac{:}));
subst=S.* (S<0);

k0=model_d.k0;
k=10.^(x(k0+1:end));
log_x=reshape(x(1:k0),model_d.nm,model_d.ncond);
back=log_x(model_d.fix_ind);
linear_x=10.^log_x;
linear_x(model_d.fix_ind)=back;

%%%%%%%%%TRANSFORM LINEAR_X%%%%%%%

meta=linear_x;
k2=repmat(k,4,1); 
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

ceq=[geq1;geq2];

