function [out_c_name,out_result_deg] = deg_processor(e_co,p_co,l_co,g_input_path)
%DEG_PROCESSOR Summary of this function goes here
%   Detailed explanation goes here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    load(strcat(g_input_path,'\process_data\ready_vector.mat'))
    disp('Data vector loaded successfully!')



    pCutoff     = p_co;
    lrCutoff    = l_co;
    existCutoff = e_co;
    
    existIndex=sum(log_data>1,2)>size(log_data,2)*existCutoff&sum(log_data_doublets>1,2)>size(log_data_doublets,2)*existCutoff;
    log_data=log_data(existIndex,:);
    log_data_doublets=log_data_doublets(existIndex,:);
    gene_name=gene_name(existIndex);
    sum(existIndex)
    
    log_data_zvalue=(log_data-repmat(mean(log_data,2),1,size(log_data,2)))./repmat(std(log_data')',1,size(log_data,2));
    log_data_doublets_zvalue=(log_data_doublets-repmat(mean(log_data_doublets,2),1,size(log_data_doublets,2)))./repmat(std(log_data_doublets')',1,size(log_data_doublets,2));
    log_data_zvalue(isnan(log_data_zvalue))=0;
    log_data_doublets_zvalue(isnan(log_data_doublets_zvalue))=0;
    
    clusterSize=max(clustering_color);
    [pvalue_total,fdr_total,logRatio_total,zvalue_total]=DEG_ranksum4cluster(clusterSize,log_data,clustering_color);
    disp('DEG ranksum Done!')
    
    
    for clusterIndex=1:clusterSize
        DEGindex(:,clusterIndex)=pvalue_total{clusterIndex}<pCutoff & logRatio_total{clusterIndex}>lrCutoff;
    end
    
    deg1 = sum(DEGindex);
    
    
    % % code to find the unique DEGs
    DEGindexOnly=zeros(size(gene_name,1),clusterSize);
    clusterOrder=1:clusterSize;
    for clusterIndex=1:size(clusterOrder,2)
        clusterIndex=clusterOrder(clusterIndex);
        for i=1:size(gene_name,1)
            DEGindexOnly(i,clusterIndex)=DEGindex(i,clusterIndex) && sum(DEGindex(i,clusterOrder),2)==1;
            if DEGindex(i,clusterIndex) && sum(DEGindex(i,clusterOrder),2)==2
                clusterTemp=clusterOrder(find(DEGindex(i,clusterOrder)));
                DEGindexOnly(i,clusterIndex)=logRatio_total{clusterIndex}(i)-logRatio_total{clusterTemp(clusterTemp~=clusterIndex)}(i)>lrCutoff;
            end
        end
    end
    
    uniqlo = sum(DEGindexOnly);
    
    DEGnumber=5;
    save_name_deg_map = strcat(g_input_path,'\visuals')
    mkdir(save_name_deg_map)
    drawHeatmap_DEG(save_name_deg_map,clustering_color,pvalue_total,pCutoff,logRatio_total,lrCutoff,log_data_zvalue,gene_name,DEGnumber,clustering_name,clustering_name_unique)
%     disp('DEG graph plotted and saved in Visuals directory!')
    
    cluster_number     = 1:length(clustering_name_unique);
    result_data_vector = vertcat(cluster_number,deg1,uniqlo);
    cluster_name       = clustering_name_unique;

    out_c_name         = cluster_name;
    out_result_deg     = result_data_vector;
    
    disp('DEG Analysis Completed!')

    msgbox("DEG graph plotted and saved in Visuals directory!","Success");

end

