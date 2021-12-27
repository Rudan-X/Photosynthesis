function[model_d]=load_experiment_data(model_d,molfold)

%% Input:
% model_d : MODEL
% molfold: mol per number of specified units

%% Output: MODEL
% model_d.exp_meta: storing experimentally measured metabolite concentrations in the specified unit

% Measured metabolite concentrations from the paper:
% Annunziata MG, Apelt F, Carillo P, Krause U, Feil R, Mengin V, et al. Getting
% back to nature: A reality check for experiments in controlled environments. Journal
% of Experimental Botany. 2017;68(16):4463{4477. doi:10.1093/jxb/erx220

% tables_meta: cell data storing the tables for the 5 experiments
% 1: fluorescence square-wave: Percival Square DLI 7 (sheet 3) 
% 2: fluorescence sinusoidal: Percival Sin DLI 7  (sheet 5)
% 3: LED square-wave: LED SQUARE 7 denser sampling  (sheet 6)
% 4: LED sinusoidal: LED Sin DLI 7 denser sampling  (sheet 7)
% 5: sunlight: Greenhouse Eq 2012  (sheet 1)


% manually copyied units from each sheet:
imp=load('../data/units_cell.mat');
units_cell=imp.units_cell;



imp=load('../data/time_inter.mat'); % time_inter: time intervals in each experiments
time_inter=imp.time_inter;

imp=load('../data/metainfo.mat'); % metainfo: the corresponding metabolites names in the tables and in the model, 
metainfo=imp.metainfo;
% with information about NAF (in compartments) and chlorophyll content


imp=load('../data/exp_data_full.mat');
tables_meta=imp.tables_meta;
full_meta = cell(5,1);
for cond = 1:5
    
    % Several measurements are taken for each time points. The average
    % value is calculated
    if cond~=5
        interv=find(~isnan(table2array(tables_meta{cond}(:,1))));
    else
        interv=find(table2array(tables_meta{cond}(:,2))~="");
    end    
    
    full_meta{cond,1} = [];
    full_meta{cond,2} = [];
    
    % find the metabolites in the tables which are present in the model:
    [~,meta_ind]=ismember(metainfo.names,string(tables_meta{cond}.Properties.VariableNames));
    
    % some of them are in micro mol/ gFW, some in nano mol/ gFW, those in nano mol is
    % divided by 1e3 in order to be converted to micro mol
    unit_vect=units_cell{cond}(2,meta_ind);
    unit_fold=ones(length(unit_vect),1);
    index = contains(unit_vect,'nmol');
    unit_fold(index)=unit_fold(index)/1000;
    
%     fprintf(' Condition: %d\n',cond);
%     check=[unit_fold,unit_vect']
    

    % for each time point, the mean and standard deviation of the measurements are taken
    for i=1:length(interv)
        
        temp1=[]; % it stores the mean of table vales according to metainfo.names within a time interval
        temp2=[];
        
        if i~=length(interv)
            range=interv(i):interv(i+1)-1;
        else
            range=interv(i):size(tables_meta{cond},1);
        end
        
        % because of different recorded format
        if cond~=5
            temp1=nanmean(table2array(tables_meta{cond}(range,meta_ind)));
            temp2=nanstd(table2array(tables_meta{cond}(range,meta_ind)),0,1);
            
        else 
            temp1=nanmean(str2double(table2array(tables_meta{cond}(range,meta_ind))));
            temp2=nanstd(str2double(table2array(tables_meta{cond}(range,meta_ind))),0,1);
        end
        
        % Change NAF to those metabolite present in both cytosol and
        % chloroplast, expanding the vector
        ct=length(temp1);
        temp1(length(metainfo.NAF))=0;
        unit_fold(length(metainfo.NAF))=0;
        temp1(ct+1:length(metainfo.NAF))=temp1([2,3,6,7,11,16]);
        unit_fold(ct+1:length(metainfo.NAF))=unit_fold([2,3,6,7,11,16]);
        
        temp2(length(metainfo.NAF))=0;
        temp2(ct+1:length(metainfo.NAF))=temp2([2,3,6,7,11,16]);
        
        
        %%%%CHANGING UNITS %%%%%%%%%%%
        % INITIAL UNITS:  micro mol/gFW  or nano mol/gFW
        % after multiplying by unit_fold: all micro mol/gFW
        
        % micro mol/gFW * 1gFW/0.8mg Chl * 1mg Chl/ 24 micro L cytosol=
        % mol/L cytosol / molfold
        % micro mol/gFW * 1gFW/0.8mg Chl * 1mg Chl/ 70 micro L plastid=
        % mol/L plastid / molfold
        
        % DEFAULT FINAL UNIT: mol/L = M
        % SPECIFIED FINAL UNIT: M / (mol per number of specified units)
        
        full_meta{cond,1}(i,:)=temp1/0.8./metainfo.chl'.*unit_fold'/molfold;
        full_meta{cond,2}(i,:)=temp2/0.8./metainfo.chl'.*unit_fold'/molfold;
 
    end
    
    % append the time points
    full_meta{cond,1}(:,end+1)=time_inter{cond,1};
    full_meta{cond,2}(:,end+1)=time_inter{cond,2};
end


%% Add the data into the model
full_meta{6,1}=metainfo.namesmodel;
model_d.exp_meta=full_meta;

% Avoid doubled time intervals :
for cond=1:5
    [U, I] = unique(model_d.exp_meta{cond,1}(:,end), 'first');
    dup = 1:length(model_d.exp_meta{cond,1}(:,end));
    dup(I) = [];
    while ~isempty(dup)
        r=[dup(1)-1,dup(1)];
        for col=1:2
            temp=model_d.exp_meta{cond,col};
            temp(r(1),:)=nanmean(model_d.exp_meta{cond,col}(r,:)); % the first one is changed by the mean vector
            temp(r(2:end),:)=[]; % the rest of duplicated ones are removed
            model_d.exp_meta{cond,col}=temp;
        end
        
        [U, I] = unique(model_d.exp_meta{cond,1}(:,end), 'first');
        dup = 1:length(model_d.exp_meta{cond,1}(:,end));
        dup(I) = [];
    end
end


