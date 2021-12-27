mat=zeros(5,23);
for cond=1:5
    ind=find(model_d.exp_meta{cond,1}(:,end)==0);
    mat(cond,:)=model_d.exp_meta{cond,1}(1,:);
end