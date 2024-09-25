close all
clear all
addpath('\\srv2\mnb\KKF\_Common\Kozgazdasagi_Modellezesi_Foosztaly\Egyeb\2024\DFM\IRIS-Toolbox-Release-20221026');
iris.startup;

%% read data

close all
clear all

for ev=2024:2024
for honap=9:9


disp('Loading data ...')
dm = databank.fromCSV('HU_monthly.csv', 'Dateformat', 'YYYY/MMM',"Delimiter",";");
dq = databank.fromCSV('HU_quarterly.csv',"Delimiter",";");

negyedev=floor((honap-1)/3)+1;    

%% cut data prior 2001 
tmp = [];
list_ = dbnames(dm);
for i = 1:numel(list_)
    tmp =[tmp dm.(list_{i})];
end

% dm = dbclip(dm, mm(2014,1):get(tmp, 'last'));
 dm = dbclip(dm, mm(2014,1):mm(ev,honap));

tmp = [];
list_ = dbnames(dq);
for i = 1:numel(list_)
    tmp =[tmp dq.(list_{i})];
end

% dq = dbclip(dq, qq(2014,1):get(tmp, 'last'));

if negyedev==1
    gdphistnegyedev=4;
    gdphistev=ev-1;
else
    gdphistnegyedev=negyedev-1;
    gdphistev=ev;
end

dq = dbclip(dq, qq(2014,1):qq(gdphistev,gdphistnegyedev));

rng_hist = qq(2014,1):qq(gdphistev,gdphistnegyedev); % historical data range for GDP
rng_nowcast = rng_hist(1):rng_hist(end)+1; % historical data range for GDP plus nowcast (1 quarter ahead)
rng_nowcast_graph = rng_hist(5):rng_hist(end)+1; % range for graphs - historical data range for GDP plus nowcast (1 quarter ahead)


%% inspect available data
% dm
% dq

%% prepare database for dfm
disp('Data transformations ...')

% seasonal adjustment
%dm = dbbatch(dm,'$1','x12(dm.$0,Inf,''mode'',''m'')','namefilter','(.*)_U','fresh',false); 
% dq = dbbatch(dq,'$1','x12(dq.$0,Inf,''mode'',''m'')','namefilter','(.*)_U','fresh',false);
%%
% log of variables and make it real
exceptions = {''};
dm = dbbatch(dm,'L_$0','100*log(dm.$0)','namelist',fieldnames(dm)-exceptions,'fresh',false);
dq = dbbatch(dq,'L_$0','100*log(dq.$0)','namelist',fieldnames(dq),'fresh',false);

dm.L_rBUX=dm.L_BUX-dm.L_CPI;
dm.L_rhitel=dm.L_hitel-dm.CPI;
dm.L_rujhitel=dm.L_ujhitel-dm.CPI;

%%
% real variables


% compute growth rates (differenciálás)
dm = dbbatch(dm,'D4L_$1','diff(dm.$0,-12)','namefilter','L_(.*)','fresh',false);
dq = dbbatch(dq,'D4L_$1','diff(dq.$0,-4)','namefilter','L_(.*)','fresh',false);

% demean data
dm_ = dbbatch(dm,'D4L_$1','dm.$0-mean(dm.$0)','namefilter','D4L_(.*)','fresh',true);
dm_.R_UNEMP= dm.unempl_rate - mean(dm.unempl_rate);
%dm_.D4_GKI_fogy_biz=diff(dm.GKI_fogy_biz,-12)-mean(diff(dm.GKI_fogy_biz,-12));
dm_.D4_GKI_fogy_biz=dm.GKI_fogy_biz-mean(dm.GKI_fogy_biz);

dq_ = dbbatch(dq,'D4L_$1','dq.$0-mean(dq.$0)','namefilter','D4L_(.*)','fresh',true);

% create dfm database
dfm_input.g_GDP_OBS  = dq_.D4L_GDP;

dfm_input.g_rujhitel_1_OBS = convert(dm_.D4L_rujhitel, 'q', 'select',1);	
dfm_input.g_rujhitel_2_OBS = convert(dm_.D4L_rujhitel, 'q', 'select',2);	
dfm_input.g_rujhitel_3_OBS = convert(dm_.D4L_rujhitel, 'q', 'select',3);

dfm_input.g_rhitel_1_OBS = convert(dm_.D4L_rhitel, 'q', 'select',1);	
dfm_input.g_rhitel_2_OBS = convert(dm_.D4L_rhitel, 'q', 'select',2);	
dfm_input.g_rhitel_3_OBS = convert(dm_.D4L_rhitel, 'q', 'select',3);

dfm_input.g_real_opg_1_OBS = convert(dm_.D4L_real_opg, 'q', 'select',1);	
dfm_input.g_real_opg_2_OBS = convert(dm_.D4L_real_opg, 'q', 'select',2);	
dfm_input.g_real_opg_3_OBS = convert(dm_.D4L_real_opg, 'q', 'select',3);

dfm_input.g_hazai_szemely_1_OBS = convert(dm_.D4L_hazai_szemely, 'q', 'select',1);	
dfm_input.g_hazai_szemely_2_OBS = convert(dm_.D4L_hazai_szemely, 'q', 'select',2);	
dfm_input.g_hazai_szemely_3_OBS = convert(dm_.D4L_hazai_szemely, 'q', 'select',3);

dfm_input.g_kulfoldi_szemely_1_OBS = convert(dm_.D4L_kulfoldi_szemely, 'q', 'select',1);	
dfm_input.g_kulfoldi_szemely_2_OBS = convert(dm_.D4L_kulfoldi_szemely, 'q', 'select',2);	
dfm_input.g_kulfoldi_szemely_3_OBS = convert(dm_.D4L_kulfoldi_szemely, 'q', 'select',3);

dfm_input.g_szemely_1_OBS = convert(dm_.D4L_szemely, 'q', 'select',1);	
dfm_input.g_szemely_2_OBS = convert(dm_.D4L_szemely, 'q', 'select',2);	
dfm_input.g_szemely_3_OBS = convert(dm_.D4L_szemely, 'q', 'select',3);

dfm_input.g_hazai_teher_1_OBS = convert(dm_.D4L_hazai_teher, 'q', 'select',1);	
dfm_input.g_hazai_teher_2_OBS = convert(dm_.D4L_hazai_teher, 'q', 'select',2);	
dfm_input.g_hazai_teher_3_OBS = convert(dm_.D4L_hazai_teher, 'q', 'select',3);

dfm_input.g_kulfoldi_teher_1_OBS = convert(dm_.D4L_kulfoldi_teher, 'q', 'select',1);	
dfm_input.g_kulfoldi_teher_2_OBS = convert(dm_.D4L_kulfoldi_teher, 'q', 'select',2);	
dfm_input.g_kulfoldi_teher_3_OBS = convert(dm_.D4L_kulfoldi_teher, 'q', 'select',3);

dfm_input.g_teher_1_OBS = convert(dm_.D4L_teher, 'q', 'select',1);	
dfm_input.g_teher_2_OBS = convert(dm_.D4L_teher, 'q', 'select',2);	
dfm_input.g_teher_3_OBS = convert(dm_.D4L_teher, 'q', 'select',3);

dfm_input.g_energia_1_OBS = convert(dm_.D4L_energia, 'q', 'select',1);	
dfm_input.g_energia_2_OBS = convert(dm_.D4L_energia, 'q', 'select',2);	
dfm_input.g_energia_3_OBS = convert(dm_.D4L_energia, 'q', 'select',3);

dfm_input.g_ingatlan_1_OBS = convert(dm_.D4L_ingatlan, 'q', 'select',1);	
dfm_input.g_ingatlan_2_OBS = convert(dm_.D4L_ingatlan, 'q', 'select',2);	
dfm_input.g_ingatlan_3_OBS = convert(dm_.D4L_ingatlan, 'q', 'select',3);

dfm_input.g_ipari_termeles_1_OBS = convert(dm_.D4L_ipari_termeles, 'q', 'select',1);	
dfm_input.g_ipari_termeles_2_OBS = convert(dm_.D4L_ipari_termeles, 'q', 'select',2);	
dfm_input.g_ipari_termeles_3_OBS = convert(dm_.D4L_ipari_termeles, 'q', 'select',3);

dfm_input.g_ipari_ertekesites_1_OBS = convert(dm_.D4L_ipari_ertekesites, 'q', 'select',1);	
dfm_input.g_ipari_ertekesites_2_OBS = convert(dm_.D4L_ipari_ertekesites, 'q', 'select',2);	
dfm_input.g_ipari_ertekesites_3_OBS = convert(dm_.D4L_ipari_ertekesites, 'q', 'select',3);

dfm_input.g_ipari_belfoldi_ertekesites_1_OBS = convert(dm_.D4L_ipari_belfoldi_ertekesites, 'q', 'select',1);	
dfm_input.g_ipari_belfoldi_ertekesites_2_OBS = convert(dm_.D4L_ipari_belfoldi_ertekesites, 'q', 'select',2);	
dfm_input.g_ipari_belfoldi_ertekesites_3_OBS = convert(dm_.D4L_ipari_belfoldi_ertekesites, 'q', 'select',3);

dfm_input.g_ipari_export_ertekesites_1_OBS = convert(dm_.D4L_ipari_export_ertekesites, 'q', 'select',1);	
dfm_input.g_ipari_export_ertekesites_2_OBS = convert(dm_.D4L_ipari_export_ertekesites, 'q', 'select',2);	
dfm_input.g_ipari_export_ertekesites_3_OBS = convert(dm_.D4L_ipari_export_ertekesites, 'q', 'select',3);

dfm_input.g_epitoipar_termeles_1_OBS = convert(dm_.D4L_epitoipar_termeles, 'q', 'select',1);	
dfm_input.g_epitoipar_termeles_2_OBS = convert(dm_.D4L_epitoipar_termeles, 'q', 'select',2);	
dfm_input.g_epitoipar_termeles_3_OBS = convert(dm_.D4L_epitoipar_termeles, 'q', 'select',3);
%
dfm_input.g_kisker_1_OBS = convert(dm_.D4L_kisker, 'q', 'select',1);	
dfm_input.g_kisker_2_OBS = convert(dm_.D4L_kisker, 'q', 'select',2);	
dfm_input.g_kisker_3_OBS = convert(dm_.D4L_kisker, 'q', 'select',3);
%
dfm_input.g_empl_1_OBS = convert(dm_.D4L_empl, 'q', 'select',1);	
dfm_input.g_empl_2_OBS = convert(dm_.D4L_empl, 'q', 'select',2);	
dfm_input.g_empl_3_OBS = convert(dm_.D4L_empl, 'q', 'select',3);
%
dfm_input.g_unempl_1_OBS = convert(dm_.D4L_unempl, 'q', 'select',1);	
dfm_input.g_unempl_2_OBS = convert(dm_.D4L_unempl, 'q', 'select',2);	
dfm_input.g_unempl_3_OBS = convert(dm_.D4L_unempl, 'q', 'select',3);
%
dfm_input.g_unempl_rate_1_OBS = convert(dm_.D4L_unempl_rate, 'q', 'select',1);	
dfm_input.g_unempl_rate_2_OBS = convert(dm_.D4L_unempl_rate, 'q', 'select',2);	
dfm_input.g_unempl_rate_3_OBS = convert(dm_.D4L_unempl_rate, 'q', 'select',3);
%
dfm_input.g_oil_huf_1_OBS = convert(dm_.D4L_oil_huf, 'q', 'select',1);	
dfm_input.g_oil_huf_2_OBS = convert(dm_.D4L_oil_huf, 'q', 'select',2);	
dfm_input.g_oil_huf_3_OBS = convert(dm_.D4L_oil_huf, 'q', 'select',3);
%
dfm_input.g_rBUX_1_OBS = convert(dm_.D4L_rBUX, 'q', 'select',1);	
dfm_input.g_rBUX_2_OBS = convert(dm_.D4L_rBUX, 'q', 'select',2);	
dfm_input.g_rBUX_3_OBS = convert(dm_.D4L_rBUX, 'q', 'select',3);
%
dfm_input.g_p_food_1_OBS = convert(dm_.D4L_p_food, 'q', 'select',1);	
dfm_input.g_p_food_2_OBS = convert(dm_.D4L_p_food, 'q', 'select',2);	
dfm_input.g_p_food_3_OBS = convert(dm_.D4L_p_food, 'q', 'select',3);

%
dfm_input.g_eu_bizalom_1_OBS = convert(dm_.D4L_EU_bizalom, 'q', 'select',1);	
dfm_input.g_eu_bizalom_2_OBS = convert(dm_.D4L_EU_bizalom, 'q', 'select',2);	
dfm_input.g_eu_bizalom_3_OBS = convert(dm_.D4L_EU_bizalom, 'q', 'select',3);

%
dfm_input.g_gki_fogy_bizalom_1_OBS = convert(dm_.D4_GKI_fogy_biz, 'q', 'select',1);	
dfm_input.g_gki_fogy_bizalom_2_OBS = convert(dm_.D4_GKI_fogy_biz, 'q', 'select',2);	
dfm_input.g_gki_fogy_bizalom_3_OBS = convert(dm_.D4_GKI_fogy_biz, 'q', 'select',3);

%return
%% plot variables
% s = get(0, 'ScreenSize');
% figure()
% 
% subplot(2,2,1)
% plot(rng_nowcast_graph, dfm_input.g_GDP_OBS, 'color', 'k', 'linewidth', 2); hold on
% plot(rng_nowcast_graph, [dfm_input.g_real_opg_1_OBS dfm_input.g_real_opg_2_OBS dfm_input.g_real_opg_3_OBS], 'linewidth', 1.5);
% legend( 'GDP', '1st month', '2nd month', '3rd month', 'Box', 'off', 'location', 'southwest')
% set(gca,'FontSize',14);
% title('Reál GDP és reál OPG, yoy')
% 
% subplot(2,2,2)
% plot(rng_nowcast_graph, dfm_input.g_GDP_OBS, 'color', 'k', 'linewidth', 2); hold on
% plot(rng_nowcast_graph, [dfm_input.g_kulfoldi_teher_1_OBS dfm_input.g_kulfoldi_teher_2_OBS dfm_input.g_kulfoldi_teher_3_OBS], 'linewidth', 1.5);
% legend( 'GDP', '1st month', '2nd month', '3rd month', 'Box', 'off', 'location', 'southwest')
% set(gca,'FontSize',14);
% title('Reál GDP és külföldi tehetforgalom, yoy')
% 
% subplot(2,2,3)
% plot(rng_nowcast_graph, dfm_input.g_GDP_OBS, 'color', 'k', 'linewidth', 2); hold on
% plot(rng_nowcast_graph, [dfm_input.g_energia_1_OBS dfm_input.g_energia_2_OBS dfm_input.g_energia_3_OBS], 'linewidth', 1.5);
% legend( 'GDP', '1st month', '2nd month', '3rd month', 'Box', 'off', 'location', 'southwest')
% set(gca,'FontSize',14);
% title('Reál GDP és energia felhasználás, yoy')
% 
% subplot(2,2,4)
% plot(rng_nowcast_graph, dfm_input.g_GDP_OBS, 'color', 'k', 'linewidth', 2); hold on
% plot(rng_nowcast_graph, [dfm_input.g_rhitel_1_OBS dfm_input.g_rhitel_2_OBS dfm_input.g_rhitel_3_OBS], 'linewidth', 1.5);
% legend( 'GDP', '1st month', '2nd month', '3rd month', 'Box', 'off', 'location', 'southwest')
% set(gca,'FontSize',14);
% title('Reál GDP és hitelszerződés, yoy')

%return
%% define model variables (without OBS suffix)
disp('Declaring the model ...')
p.list={"g_GDP",...
"g_ipari_termeles_1", "g_ipari_termeles_2", "g_ipari_termeles_3", ...
"g_ipari_ertekesites_1", "g_ipari_ertekesites_2", "g_ipari_ertekesites_3", ...
"g_unempl_1", "g_unempl_2", "g_unempl_3", ...
"g_empl_1", "g_empl_2", "g_empl_3", ...
"g_gki_fogy_bizalom_1", "g_gki_fogy_bizalom_2", "g_gki_fogy_bizalom_3", ...
"g_eu_bizalom_1", "g_eu_bizalom_2", "g_eu_bizalom_3", ...
"g_szemely_1", "g_szemely_2", "g_szemely_3", ...
"g_real_opg_1", "g_real_opg_2", "g_real_opg_3", ...
"g_teher_1", "g_teher_2", "g_teher_3", ...
"g_hazai_teher_1", "g_hazai_teher_2", "g_hazai_teher_3", ...
"g_energia_1", "g_energia_2", "g_energia_3"};

% "g_hazai_teher_1", "g_hazai_teher_2", "g_hazai_teher_3", ...
% "g_ipari_belfoldi_ertekesites_1", "g_ipari_belfoldi_ertekesites_2", "g_ipari_belfoldi_ertekesites_3", ...
% "g_ipari_export_ertekesites_1", "g_ipari_export_ertekesites_2", "g_ipari_export_ertekesites_3", ...
% "g_rhitel_1", "g_rhitel_2", "g_rhitel_3", ...
% "g_epitoipar_termeles_1", "g_epitoipar_termeles_2", "g_epitoipar_termeles_3", ...
% "g_kisker_1", "g_kisker_2", "g_kisker_3", ...




  




%"g_teher_1", "g_teher_2", "g_teher_3", ...
%"g_rhitel_1", "g_rhitel_2", "g_rhitel_3", ...
%"g_gki_fogy_bizalom_1", "g_gki_fogy_bizalom_2", "g_gki_fogy_bizalom_3", ...

%"g_rBUX_1", "g_rBUX_2", "g_rBUX_3", ...
%"g_oil_huf_1", "g_oil_huf_2", "g_oil_huf_3", ...
%"g_p_food_1", "g_p_food_2", "g_p_food_3", ...
%"g_unempl_rate_1", "g_unempl_rate_2", "g_unempl_rate_3", ...
%"g_rhitel_1", "g_rhitel_2", "g_rhitel_3", ...
%"g_teher_1", "g_teher_2", "g_teher_3", ...


%
% "g_unempl_rate_1", "g_unempl_rate_2", "g_unempl_rate_3", ...


%"g_szemely_1", "g_szemely_2", "g_szemely_3", ...
%"g_ingatlan_1", "g_ingatlan_2", "g_ingatlan_3", ...
%"g_ferihegy_1", "g_ferihegy_2", "g_ferihegy_3", ...
%"g_Gtrends_munkanelk_1", "g_Gtrends_munkanelk_2", "g_Gtrends_munkanelk_3", ...
%"g_Gtrends_allasker_1", "g_Gtrends_allasker_2", "g_Gtrends_allasker_3"

% "g_rujhitel_1", "g_rujhitel_2", "g_rujhitel_3", ...
% "g_hazai_szemely_1", "g_hazai_szemely_2", "g_hazai_szemely_3", ...
% "g_kulfoldi_szemely_1", "g_kulfoldi_szemely_2", "g_kulfoldi_szemely_3", ...
% "g_kulfoldi_teher_1", "g_kulfoldi_teher_2", "g_kulfoldi_teher_3", ...



%%
m = Model.fromFile('dfm.model', 'linear', true, 'assign', p, 'AutoDeclareParameters', true);

% initial parametrization of the model -- starting point of estimation
parList = get(m, 'pList');
for i = 1:numel(parList)
   par.(parList{i}) = 0.5; % initial value of parameters for maximum likelihood estimate
end
%%
m = assign(m, par);

% inspect the model structure
get(m, 'equations')

% return
%% estimation of parameters
disp('Estimating model parameters ...')
[summary, poster, proposalCov, hess, mEst] ...
          = estimate(m, dfm_input, rng_hist, get(m, 'parameters'), 'NoSolution', 'Penalty', 'MaxFunEvals', 2e4);

summary(:, ["PosterMode", "PosterStd", "Bounds", "Info", "Start"])
%return
%% filtration including forecast
disp('Fitration of the data and nowcasting ...')

f0 = kalmanFilter(mEst, dfm_input, rng_nowcast);


% figure
% plot(rng_nowcast, [f0.Mean.g_GDP f0.Mean.g_GDP-f0.Mean.SHK_RES_g_GDP], 'linewidth', 2, 'DateTick', rng_nowcast(1):8:rng_nowcast(end));
% set(gca, 'Box', 'off')
% title('Real GDP growth (de-meaned), % yoy')
% legend('Actual', 'Predicted', 'Box', 'off')
% grid
    
%return
%% store results


%%
% load out
out = dq*{'GDP', 'L_GDP'};
out = dbmerge(out, f0.Mean*{'S', 'RES_g_GDP'});
seged=f0.Mean.g_GDP-f0.Mean.SHK_RES_g_GDP+mean(dq.D4L_GDP);
out.nowcast = seged;
out.Mean_D4L_GDP = mean(dq.D4L_GDP);
out.lambda_g_GDP=summary{1,1};
seged='dfm_outputs_'+ string(ev)+'_'+string(honap)+".csv";
databank.toCSV(out, seged);

end
end
%%
f0.Mean.S(end)*summary{1,1}+mean(dq.D4L_GDP)+f0.Mean.RES_g_GDP(end)
f0.Mean.g_GDP(end)-f0.Mean.SHK_RES_g_GDP(end)+mean(dq.D4L_GDP)
f0.Mean.g_GDP-f0.Mean.SHK_RES_g_GDP+mean(dq.D4L_GDP)
