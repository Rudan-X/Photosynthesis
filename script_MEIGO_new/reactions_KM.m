% % 1 <=> 1
% reac11=['tpic', 'rpi','rpe'];
% 
% kmf=[nan,]
% kmp=[1.32,]
% kcatf=[3406, 
% kcatb=[


% 2 <=> 2
reac22={'pgaks', 'transk', 'transk2',...
    'pgakc', 'gapdhc','fbpasec', 'ndpk','pfk2','ugpase','sps'};

kma=[0.39,0.1,0.46,0.39,0.004,0.7,1.66,0.55,0.14, 2.4];
kmb=[0.24,0,1,0.072,0.24,nan,12,0.28,5,0.1,0.8];
kmp=[0.23,0.1,0.1,0.23,nan,0.033,16,0.16,0.11,0.4];
kmq=[nan,nan,nan,nan,nan,nan,0.042,0.021,0.12,0.7];

bounds=cell(length(reac22),1);

options=optimoptions('fmincon','MaxFunctionEvaluations', 6000);
for i=6:length(reac22)
    display(i)
    reac=reac22{i};
    ind=find(contains(model_d.rxns,reac));
    if i==2
        ind=ind([1:5,11:15]);
    elseif i==3
        ind=ind([2:6,8:12]);
    end
    kms=[kma(i),kmb(i),kmp(i),kmq(i)];
    nonlcon = @(x)mycon(x,kms);
    bounds{i}=zeros(length(ind),2);
    for k=1:length(ind)
        display(k)
        x0=[];A=[]; b=[]; Aeq=[]; beq=[];
        lb=1e-6*ones(length(ind),1);
        ub=1e6*ones(length(ind),1);
        optmin = fmincon(@(x)x(k),rand(10,1),A,b,Aeq,beq,lb,ub,nonlcon);
        bounds{i}(k,1)=optmin(k);
        optmax = fmincon(@(x)(-x(k)),rand(10,1),A,b,Aeq,beq,lb,ub,nonlcon);
        bounds{i}(k,2)=optmax(k);
    end
end



