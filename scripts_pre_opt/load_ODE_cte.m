function[model_d,metacte]=load_ODE_cte(model_d,molfold,timefold)
%% Input:
% model_d : MODEL
% molfold: mol per number of specified units

%% Output:
% cbcmeta: structure storing fiexed metabolite concentrations in the specified
% units and their indices in the model

% model_d.irradiance: structure storing the fixed flux for the reaction
% which import the photons

%%
model_d.constant=zeros(length(model_d.mets),2);

metacte=struct();

% list the metabolites with known information
metacte.names=string({'Pi[s]','Pi[c]','CO2[s]', 'O2[s]','NADp','NADH'});
% The amount of the corresponding compounds given in miliM, then converted
% to M
metacte.content=zeros(length(metacte.names),2);

% day concentration
metacte.content(:,1)=[2; 31;  0.33; 0.245; 0.6; 0.005]*(1e-3);

% dark concentration
metacte.content(:,2)=[5; 17;  0.33; 0.245; 1;   0.007]*(1e-3);


%%%%CHANGING UNITS %%%%%%%%%%%
% DEFAULT FINAL UNIT: mol/L = M
% SPECIFIED FINAL UNIT: M / (mol per number of specified units)

metacte.content=metacte.content/molfold;

[~,metacte.ind]=ismember(metacte.names,string(model_d.mets));
model_d.constant(metacte.ind,:)=metacte.content;


%% Irradiance patterns of different conditions
% Stored as the same format as model_d.exp_meta for every time point

model_d.irradiance=cell(5,1);

% condition 1: blue square symbols (fluorescence square):
model_d.irradiance{1}=[[0;repmat(165,11,1);zeros(11,1)],model_d.exp_meta{1,1}(:,end)];

% condition 2: blue triangle symbols (fluorescence sinusoidal):
model_d.irradiance{2}=[[0;0;50;75;100;140;200;240;200;120;75;30;0;zeros(11,1)],model_d.exp_meta{2,1}(:,end)];

% condition 3: gray square symbols (LED square):
model_d.irradiance{3}=[[0;repmat(160,10,1);zeros(8,1)],model_d.exp_meta{3,1}(:,end)];

% condition 4: gray triangle symbols (LED sinusoidal):
model_d.irradiance{4}=[[0;50;75;150;210;250;200;120;60;20;0;zeros(8,1)],model_d.exp_meta{4,1}(:,end)];

% condition 5: orange circle symbols (natural light):
model_d.irradiance{5}=[[0;60;85;162.5;100;225;205;160;163;230;162;230;60;100;30;0;zeros(14,1)],model_d.exp_meta{5,1}(:,end)];

%changing units
% light input [lumen] ,flux given in units: MICROmol/s/m2

% micromol/s/m2 * 1m2 leaf /(0.5*1e-3*58.5%*9.5%)m^3 stroma * 1m^3/1000L 
% = 10e-6 mol/L/s
% units: mol/s/L [lumen]
for cond=1:5
    model_d.irradiance{cond}(:,1)=areatoliter(model_d.irradiance{cond}(:,1))* 1e-6 * timefold / molfold;
end

