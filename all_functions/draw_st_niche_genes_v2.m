  function draw_st_niche_genes(current_path,spot_size)
%     current_path = 'C:\Users\ja_ha\OneDrive\Documents\MATLAB\Chapter_5\PDAC_A\ST_1\Results_05';
    comb_list_dir = strcat(current_path,'\data\deglist');
    load_filename = strcat(current_path,'\process_data\ready_vector.mat')
    load(load_filename)
    
    Files = struct2cell(dir(fullfile(comb_list_dir, '*.mat')))';
    
    for i = 1:length(Files(:,1))
        comb_files(i,1) = convertCharsToStrings(strcat(Files(i,2),'\',Files(i,1)));
    end
    
    total_combination = convertCharsToStrings(Files(:,1));
    total_combination = extractBefore (total_combination,'.mat');
    
    %%
    % load('n0')
    load(strcat(current_path,'\data\SVM_bestMatch.mat'))
    
    col_names = st_cordinates';
    
    cordinates = split(col_names,"x");
    x=str2double(cordinates(:,1))';
    y=str2double(cordinates(:,2))';
    x2 = {};
    y2 = {};
    
    for i=1:length(total_combination(:,1))
        ind = find(strcmp(artificialDoubletsCombiUnique,total_combination(i)));
        imp_ind = find(bestMatch==ind);
        col_data = cordinates(imp_ind',:);
    %     cordinates2 = split(col_data,"x");
        x2{i,1} = str2double(col_data(:,1))';
        y2{i,1} = str2double(col_data(:,2))';
        clear ind imp_ind col_data cordinates2
    end
    %% Phase 2
    % this part plot all the indexes of genes
    
    combinationss = length(total_combination(:,1));
            
    existCutoff = 0.01;
    existIndex=sum(log_data>1,2)>size(log_data,2)*existCutoff&sum(log_data_doublets>1,2)>size(log_data_doublets,2)*existCutoff;
    log_data=log_data(existIndex,:);
    log_data_doublets=log_data_doublets(existIndex,:);
    gene_name=gene_name(existIndex);
    sum(existIndex)
    
    log_data_zvalue=(log_data-repmat(mean(log_data,2),1,size(log_data,2)))./repmat(std(log_data')',1,size(log_data,2));
    log_data_doublets_zvalue=(log_data_doublets-repmat(mean(log_data_doublets,2),1,size(log_data_doublets,2)))./repmat(std(log_data_doublets')',1,size(log_data_doublets,2));
    
    log_data_zvalue(isnan(log_data_zvalue))=0;
    log_data_doublets_zvalue(isnan(log_data_doublets_zvalue))=0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    comb_files
    for j=1:combinationss
    
    %     niche_genes = load(comb_files(j,1));
        load(comb_files(j,1))
        total_niche_genes{j,1} = niche_geness;

        counter = 1;
        genes_to_draw = 10;

        if length(niche_geness) < 10
            genes_to_draw = length(niche_geness);
        end
        
        if length(niche_geness) > 0

            for i=1:genes_to_draw
                current_gene = niche_geness(i,1);
                temp_index = find(strcmp(gene_name,current_gene));
                if ~isempty(temp_index)
                    index(counter,1) = temp_index;
                    counter = counter + 1;
                end
            end
        
            c_mat{j,1} = log_data_doublets_zvalue(index,:);
        
            cordinates3 = split(col_names,"x");
            x3 = str2double(cordinates(:,1))';
            y3 = str2double(cordinates(:,2))';

        end
        clear index niche_geness temp_index index_data index_data_max
        clear temp_col_names_2_2
    end
    
    % clearvars -except t_name
    
    %% Ploting section
    clearvars -except total_niche_genes current_path t_name x x2 x3 y y2 y3 c_mat total_combination niche_genes spot_size
    sz=spot_size;
    % for i=1:1
    for i=1:length(total_combination)

        niche_geness = total_niche_genes{i,1}
        c_mat2 = cell2mat(c_mat(i,:));
        fig1 = figure;
        scatter(x,y,sz,'filled','s')
        hold on
        scatter(cell2mat(x2(i)),cell2mat(y2(i)),sz,'filled','sr')

        view(90,90)
        x0=0;
        y0=0;
        width  = 518;
        height = 518;
        set(gcf,'position',[x0,y0,width,height])
        xlim([0 max(x)])
        ylim([0 max(y)])

        t = strcat("Predicted spots for ",total_combination(i));
        title(t)
    %     title(total_combination(i))
        mkdir(strcat(current_path,'\ST_visual\Neibour_specific\',total_combination(i)))
        fig1_name = strcat(current_path,'\ST_visual\Neibour_specific\',total_combination(i),'\best_match_',total_combination(i),'_',num2str(i),'.png')
        saveas(fig1,fig1_name)
        close(fig1)
        
        genes_to_draw = 10;

        if length(niche_geness) < 10
            genes_to_draw = length(niche_geness);
        end

        for j=1:genes_to_draw
            fig2 = figure;
            scatter(x,y,sz,'filled','s')
            hold on    
            scatter(x3,y3,sz,c_mat2(j,:),'filled','s')
            colormap(gca,'Jet')
            colorbar

            view(90,90)
            x0=0;
            y0=0;
            width  = 510;
            height = 450;
            set(gcf,'position',[x0,y0,width,height])
            xlim([0 max(x)])
            ylim([0 max(y)])

            t = strcat({'Combination: '},total_combination(i),{' & Gene: '},niche_geness(j));
            title(t)  
            fig2_name = strcat(current_path,'\ST_visual\Neibour_specific\',total_combination(i),'\',num2str(j),'_niche_',total_combination(i),'_',niche_geness(j),'.png')
            saveas(fig2,fig2_name)
            close(fig2)
        end
    end

    msgbox("Neighbor specific genes graph plotted!","Success");
end