function[model_d,model_orig]=decomp_model(model,enzymes,cat)
%having the Stoichiometric matrix (S), the vector with catalytic info (cat)
%we can first get the dimension of the expanded stoichiometric matrix S_exp
%m metabolite and n reactions in the original S. p reactions are catalyzed
%by enzymes

model_orig=model;

%%
% Add the additional enzymes to the metabolite list
for e=1:length(enzymes.name)
    if sum(contains(model.mets,enzymes.name{e}))==0
        model.mets{length(model.mets)+1}=enzymes.name{e};
    end
end


% Divide the set of reactions into: catalyzed by enzymes vs not catalyzed
ind_simp=[];
ind_decomp=[];
% cat_new=cell(0);
for r=1:length(model.rxns)
    if isempty(cat{r})
        ind_simp=[ind_simp;r];
    else
        ind_decomp=[ind_decomp;r];
%         cat_new{length(cat_new)+1}=cat{r};
    end
end



% the dimension of the expanded matrix should be:

mat1=model_orig.S(:,ind_simp);

new_n=0;
new_m=0;
for j=1:length(model_orig.rxns)
    if ~isempty(cat{j})
        ind_s=find(model_orig.S(:,j)<0);
        ind_p=find(model_orig.S(:,j)>0);
        new_n=new_n+1+length(ind_s)+length(ind_p);
        new_m=new_m+length(ind_s)+length(ind_p);
    end
end



% mE=max(cat);% number of enzymes
% new_m=mI+mE+new_m;


enzy_ind=[]; %for each reaction, the enzymes(index) catalyzing it is recorded
% meta_enzy=zeros(new_m,1); % for each enzymes complex form, the enzymes(index) it belongs to, pure metabolite will have value 0
cxe=zeros(size(model.S,1),length(enzymes.name));

S_addit=zeros(0);
current_I=length(model.mets);
current_J=0;
rev_decomp=ones(new_n,1);
exp_rxns={};
%%
for j=1:length(model.rxns)
    if ~isempty(cat{j})
        ind_s=find(model.S(:,j)<0);
        ind_p=find(model.S(:,j)>0);

        [~,ind_enz]=ismember(cat{j},string(model.mets));
        [~,ind_temp]=ismember(cat{j},enzymes.name);
        ct=0;
        cxe(ind_enz,ind_temp)=1;
        for s=1:length(ind_s)
            ct=ct+1;
            %new column for the new elementary reaction
            current_J=current_J+1;
            exp_rxns(current_J)={append(model.rxns{j},num2str(ct))};
            enzy_ind(current_J)=ind_temp;

 
            %when the reactant is a substrate:
            %S + E -> SE
            if s==1
                S_addit(ind_enz,current_J)=-1;
                temp=strcat(model.mets{ind_enz},"_",model.mets{ind_s(s)});
                [~,ind]=ismember(temp,string(model.mets));
                if ind==0
                     %new row for the new reactant enzymes complex if it
                     %does not exist before
                    current_I=current_I+1;
                    ind_complex=current_I;
                else
                    % if it exist, then tthe index is assigned
                    ind_complex=ind;
                end
            else
                S_addit(ind_complex,current_J)=-1; % from previous round
                temp=strcat(string(model.mets(ind_complex)),"_",model.mets{ind_s(s)});
                [~,ind]=ismember(temp,string(model.mets));
                if ind==0
                    current_I=current_I+1;
                    ind_complex=current_I;
                else
                    ind_complex=ind;
                end
            end
            model.mets{ind_complex}=temp;
            S_addit(ind_s(s),current_J)=-1;
            S_addit(ind_complex,current_J)=1;
            cxe(ind_complex,ind_temp)=1;
        end
        
        ct=ct+1;
        %new column for the new elementary reaction: SE ->PE
        current_J=current_J+1;
        exp_rxns(current_J)={append(model.rxns{j},num2str(ct))};
        enzy_ind(current_J)=ind_temp;

        S_addit(current_I,current_J)=-1;
        S_addit(current_I+1,current_J)=1;
        
        allprod=strjoin(string(model.mets(ind_p)),'_');
        model.mets{current_I+1}=strcat(model.mets{ind_enz},"_",allprod);
        cxe(current_I+1,ind_temp)=1;
        
        for p=1:length(ind_p)
            ct=ct+1;
            %new column for the new elementary reaction
            current_J=current_J+1;
            exp_rxns(current_J)={append(model.rxns{j},num2str(ct))};
            enzy_ind(current_J)=ind_temp;
            %enzy_ind(current_J)=cat(j);
            %the product is not reversible
            if model.rev(j)==0
                rev_decomp(current_J)=0;
            end
            %new row for the new product enzymes complex
            current_I=current_I+1;
            cxe(current_I,ind_temp)=1;
            %%when the reactant is a product:
            %PE -> P + E
            if p==length(ind_p)
                S_addit(ind_enz,current_J)=1;
            else
                S_addit(current_I+1,current_J)=1;
                produ=strcat('_',model.mets{ind_p(p)});
                model.mets{current_I+1}=erase(string(model.mets(current_I)),produ);
            end
            S_addit(ind_p(p),current_J)=1;
            S_addit(current_I,current_J)=-1;
        end
    end
end
%%
mat1(size(S_addit,1),end)=0; %expand the metabolite
final_stoich=[mat1,S_addit];


% Building the structure
model_d=struct();
model_d.S=final_stoich;
model_d.rxns=[model.rxns(ind_simp);exp_rxns'];
model_d.mets=model.mets;
model_d.c=zeros(size(final_stoich,2),1);
model_d.ub=100000*ones(size(final_stoich,2),1);
model_d.lb=zeros(size(final_stoich,2),1);
model_d.rev=[model.rev(ind_simp);rev_decomp];
model_d.lb(model_d.rev==1)=-10000;
model_d.b=zeros(size(final_stoich,1),1);
model_d.cata=[zeros(length(ind_simp),1);enzy_ind'];
model_d.cata(find(contains(string(model_d.rxns),'_b')))=-model_d.cata(find(contains(string(model_d.rxns),'_b')));
model_d.enzyme_names=[string(enzymes.name),string(enzymes.EC)];
model_d.cxe=cxe;

% include the enzymes which are not catalytic (structural or inactive form)

model_d.enzyme_names(end+1,1)='photos2';
[~,metsind]=ismember(string(enzymes.noncata(1:20)),string(model_d.mets));
model_d.cxe(metsind,end+1)=1;

model_d.enzyme_names(end+1,1)='cytochrome';
[~,metsind]=ismember(string(enzymes.noncata(21:31)),string(model_d.mets));
model_d.cxe(metsind,end+1)=1;

model_d.enzyme_names(end+1,1)='photos1';
[~,metsind]=ismember(string(enzymes.noncata(32:36)),string(model_d.mets));
model_d.cxe(metsind,end+1)=1;

model_d.enzyme_names(end+1,1)='chlorophyll';
[~,metsind]=ismember(string(enzymes.noncata(63:end)),string(model_d.mets));
model_d.cxe(metsind,end+1)=1;


[~,metsind]=ismember(string(enzymes.noncata(37)),string(model_d.mets));
[~,enzind]=ismember('atpsac',string(enzymes.name));
model_d.cxe(metsind,enzind)=1;

model_d.enzyme_names(end+1,1)='rubisco_activase';
[~,metsind]=ismember(string(enzymes.noncata(38:39)),string(model_d.mets));
model_d.cxe(metsind,end+1)=1;

model_d.enzyme_names(end+1,1)='rubisco';
[~,metsind]=ismember(string(enzymes.noncata(40:57)),string(model_d.mets));
model_d.cxe(metsind,end+1)=1;

[~,metsind]=ismember(string(enzymes.noncata(58:59)),string(model_d.mets));
[~,enzind]=ismember('psbsh6zx',string(enzymes.name));
model_d.cxe(metsind,enzind)=1;

[~,metsind]=ismember(string(enzymes.noncata(60)),string(model_d.mets));
[~,enzind]=ismember('vdeh3',string(enzymes.name));
model_d.cxe(metsind,enzind)=1;

[~,metsind]=ismember(string(enzymes.noncata(61)),string(model_d.mets));
[~,enzind]=ismember('stn7a',string(enzymes.name));
model_d.cxe(metsind,enzind)=1;

[~,metsind]=ismember(string(enzymes.noncata(62)),string(model_d.mets));
[~,enzind]=ismember('pph1a',string(enzymes.name));
model_d.cxe(metsind,enzind)=1;



