clear all
% mkdir Empirical_Result/Graph
load ('Empirical_Result/macrofactor/IPCADATA_FNW36_RNKDMN_CON_K5_macrosLASSOFLEX.mat');
%%
characteristics_names={'beta','e2p','bm','q','ltrev','a2me','strev','mktcap','cto','oa','roe','mom','turn','noa','w52h','idiovol','ato','dpi2a','invest','pm','rna','suv','roa','freecf','ol','c','intmom','bidask','assets','lev','prof','s2p','sga2m','d2a','fc2y','pcm','constant'}

%% Plot barchart
GBGD = [GBGD_GB,GBGD_GD] 

for i=2:2
b=bar(GBGD(:,i:i));
xticks(1:37)
ylim([-0.8,0.8]);
yticks(-0.8:0.2:0.8)
 title(['Factor ' num2str(i)],'FontSize',16,'FontName', 'Times New Roman');
set(b, 'FaceColor', 'k'); 
xticklabels({'beta','e2p','bm','q','ltrev','a2me','strev','mktcap','cto','oa','roe','mom','turn','noa','w52h','idiovol','ato','dpi2a','invest','pm','rna','suv','roa','freecf','ol','c','intmom','bidask','assets','lev','prof','s2p','sga2m','d2a','fc2y','pcm','constant'});
gca.XAxis.FontSize = 4
xtickangle(90);

saveas(b, sprintf('Empirical_Result/Graph/''macro_FIG%d.png',i));
end
%%