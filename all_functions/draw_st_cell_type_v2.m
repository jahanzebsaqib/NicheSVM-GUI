function draw_st_cell_type(current_path,spot_size)

    load_filename = strcat(current_path,'\process_data\ready_vector.mat')
    load(load_filename)
    

    load_filename = strcat(current_path,'\data\use_parameters.mat')
    load(load_filename)


    cell_names = clustering_name_unique';

    existCutoff = n_e_co; pCutoff = n_p_co; lrCutoff = n_l_co;

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

    for clusterIndex=1:clusterSize
        DEGindex(:,clusterIndex)=pvalue_total{clusterIndex}<pCutoff & logRatio_total{clusterIndex}>lrCutoff;
    end
    deg1 = sum(DEGindex)
    total_deg = sum(deg1)
    
    for total_cells=1:length(cell_names)
        markers_genes = gene_name(find(DEGindex(:,total_cells)));
    
        for genes_loops=1:length(markers_genes)
            z_data(genes_loops,:) = log_data_doublets_zvalue(find(strcmp(gene_name,markers_genes(genes_loops,:))),:);
        end
        z_mean(total_cells,:) = mean(z_data,1);
    end    
        
    clearvars -except z_mean curr_path cell_names st_cordinates current_path spot_size
    
    %%
    
    col_names = st_cordinates';
    
    cordinates = split(col_names,"x");
    x=str2double(cordinates(:,1))';
    y=str2double(cordinates(:,2))';
    sz=spot_size;
    
%     load_filename = strcat(current_path,'\process_data\ready_vector.mat')
%     load(load_filename , 'clustering_name_unique')
    

    mkdir(strcat(current_path,'\ST_visual\Cell_type'))
            
%     mkdir 'figs'

    for plot_loop=1:length(cell_names)
        fig2 = figure;
        scatter(x,y,sz,z_mean(plot_loop,:),'filled','s')
        colormap(gca,'Jet')
        colorbar
        title(cell_names(plot_loop))
        view(90,90)
        x0=0;
        y0=0;
        width  = 510;
        height = 450;
        set(gcf,'position',[x0,y0,width,height])
        xlim([0 max(x)])
        ylim([0 max(y)])
        fig2_name = strcat(current_path,'\ST_visual\Cell_type\',num2str(plot_loop),'_',cell_names(plot_loop),'.png')
        saveas(fig2,fig2_name)
        close(fig2)
    end

    msgbox("Cell type expression graph plotted!","Success");

end

