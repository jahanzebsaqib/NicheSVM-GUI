function drawHeatmap_DEG(outputFolder,clustering_color,pvalue_total,pCutoff,logRatio_total,lrCutoff,log_data_zvalue,gene_name,DEGnumber,clustering_name,clustering_name_unique)

load 'colormap_2to19grey.mat'

clusterSize=max(clustering_color);
DEGindex=zeros(size(gene_name,1),clusterSize);
for clusterIndex=1:clusterSize
    DEGindex(:,clusterIndex)=pvalue_total{clusterIndex}<pCutoff & logRatio_total{clusterIndex}>lrCutoff;
end

DEGindexOnly=zeros(size(gene_name,1),clusterSize);
clusterOrder=1:clusterSize;
for clusterIndex=1:size(clusterOrder,2)
    for i=1:size(gene_name,1)
        DEGindexOnly(i,clusterIndex)=DEGindex(i,clusterIndex) && sum(DEGindex(i,clusterOrder),2)==1;
        if DEGindex(i,clusterIndex) && sum(DEGindex(i,clusterOrder),2)==2
            clusterTemp=clusterOrder(find(DEGindex(i,clusterOrder)));
            DEGindexOnly(i,clusterIndex)=logRatio_total{clusterIndex}(i)-logRatio_total{clusterTemp(clusterTemp~=clusterIndex)}(i)>lrCutoff;
        end
    end
end

geneIndex=[];geneIndex2=[];cellIndex=[];
for clusterIndex=1:clusterSize
    geneIndexTemp=find(DEGindex(:,clusterIndex));
    [~,sortIndex]=sort(logRatio_total{clusterIndex}(geneIndexTemp),'descend');
    geneIndexTemp=geneIndexTemp(sortIndex);
    for ii=1:size(geneIndexTemp,1)
        if sum(geneIndex==geneIndexTemp(ii))==0
            geneIndex=[geneIndex;geneIndexTemp(ii)];
        end
    end
    if size(geneIndexTemp,1)>DEGnumber
        geneIndexTemp2=geneIndexTemp(1:DEGnumber);
    else
        geneIndexTemp2=geneIndexTemp;
    end
    geneIndex2=[geneIndex2;geneIndexTemp2];
    cellIndexTemp=find(clustering_color==clusterIndex);
    [~,sortIndex]=sort(mean(log_data_zvalue(geneIndexTemp,cellIndexTemp),1),'descend');
    cellIndexTemp=cellIndexTemp(sortIndex);
    cellIndex=[cellIndex cellIndexTemp];
end
close all
figure(1)
ax(1)=subplot(2,1,1);
imagesc(log_data_zvalue(geneIndex2,cellIndex))
xticks([])
yticks([1:size(geneIndex2,1)])
yticklabels(gene_name(geneIndex2))
caxis([-3 3])
colormap jet
set(gca, 'Fontsize', 7)
set(gcf, 'Position', [100, 100, 400, 600])
ax(2)=subplot(2,1,2);
imagesc(clustering_color(cellIndex))
% xticks([])
cluster_portion = clustering_color(cellIndex);
[unique_vals,unique_vals_index] = unique(cluster_portion);
for ini_loop=1:length(unique_vals_index)-1
    name_indexes(ini_loop) = round((unique_vals_index(ini_loop) + unique_vals_index(ini_loop+1))/2);
end
name_indexes(ini_loop+1) = floor((unique_vals_index(ini_loop+1) + length(cellIndex))/2);
name_indexes = name_indexes-1;
cluster_disp_name = clustering_name_unique(clusterOrder);
xticks([name_indexes])
xticklabels(cluster_disp_name)
xtickangle(45)

yticks([])
caxis([0 clusterSize])
colormap(ax(2),mycmap2to19grey{clusterSize-1})
subplot(2,1,2)
p = get(gca, 'Position');
p(4) = p(4) / 8;
set(gca, 'Position', p);
subplot(2,1,1)
p = get(gca, 'Position');
p_diff = p(4) * 1.2;
p(4) = p(4) + p_diff;
p(2) = p(2) - p_diff;
set(gca, 'Position', p);
saveas(gcf,[outputFolder,'/heatmap_cluster_pvalue',num2str(pCutoff),'lr',num2str(lrCutoff),'top',num2str(DEGnumber),'.pdf'])
saveas(gcf,[outputFolder,'/heatmap_cluster_pvalue',num2str(pCutoff),'lr',num2str(lrCutoff),'top',num2str(DEGnumber),'.fig'])






geneIndex=[];geneIndex2=[];cellIndex=[];
clusterOrder=find(sum(DEGindexOnly)>0);
for clusterIndex=1:size(clusterOrder,2)
    clusterIndex=clusterOrder(clusterIndex);
    geneIndexTemp=find(DEGindexOnly(:,clusterIndex));
    [~,sortIndex]=sort(logRatio_total{clusterIndex}(geneIndexTemp),'descend');
    geneIndexTemp=geneIndexTemp(sortIndex);
    geneIndex=[geneIndex;geneIndexTemp];
    if size(geneIndexTemp,1)>DEGnumber
        geneIndexTemp2=geneIndexTemp(1:DEGnumber);
    else
        geneIndexTemp2=geneIndexTemp;
    end
    geneIndex2=[geneIndex2;geneIndexTemp2];
    cellIndexTemp=find(clustering_color==clusterIndex);
    [~,sortIndex]=sort(mean(log_data_zvalue(geneIndexTemp,cellIndexTemp),1),'descend');
    cellIndexTemp=cellIndexTemp(sortIndex);
    cellIndex=[cellIndex cellIndexTemp];
end
close all
figure(1)
ax(1)=subplot(2,1,1);
imagesc(log_data_zvalue(geneIndex2,cellIndex))
xticks([])
yticks([])
yticks([1:size(geneIndex2,1)])
yticklabels(gene_name(geneIndex2))
caxis([-3 3])
colormap jet
set(gca, 'Fontsize', 7)
set(gcf, 'Position', [100, 100, 400, 600])
ax(2)=subplot(2,1,2);
imagesc(clustering_color(cellIndex))
cluster_portion = clustering_color(cellIndex);
[unique_vals,unique_vals_index] = unique(cluster_portion);
for ini_loop=1:length(unique_vals_index)-1
    name_indexes(ini_loop) = round((unique_vals_index(ini_loop) + unique_vals_index(ini_loop+1))/2);
end
name_indexes(ini_loop+1) = floor((unique_vals_index(ini_loop+1) + length(cellIndex))/2);
name_indexes = name_indexes-1
cluster_disp_name = clustering_name_unique(clusterOrder);
xticks([name_indexes])
xticklabels(cluster_disp_name)
xtickangle(45)

yticks([])
caxis([0 clusterSize])
colormap(ax(2),mycmap2to19grey{clusterSize-1})
subplot(2,1,2)
p = get(gca, 'Position');
p(4) = p(4) / 8;
set(gca, 'Position', p);
subplot(2,1,1)
p = get(gca, 'Position');
p_diff = p(4) * 1.2;
p(4) = p(4) + p_diff;
p(2) = p(2) - p_diff;
set(gca, 'Position', p);
saveas(gcf,[outputFolder,'/heatmap_clusterOnly_pvalue',num2str(pCutoff),'lr',num2str(lrCutoff),'top',num2str(DEGnumber),'.pdf'])
saveas(gcf,[outputFolder,'/heatmap_clusterOnly_pvalue',num2str(pCutoff),'lr',num2str(lrCutoff),'top',num2str(DEGnumber),'.fig'])

