function nicheSVM_processor(n_e_co , n_p_co , n_l_co , n_clusters , n_g_input_path)
%NICHESVM_PROCESSOR Summary of this function goes here
%   Detailed explanation goes here

    try
        load(strcat(n_g_input_path,'\process_data\ready_vector.mat'))
        disp('Data vector loaded successfully!')
        existCutoff = n_e_co;
        existIndex=sum(log_data>1,2)>size(log_data,2)*existCutoff&sum(log_data_doublets>1,2)>size(log_data_doublets,2)*existCutoff;
        log_data=log_data(existIndex,:);
        log_data_doublets=log_data_doublets(existIndex,:);
        gene_name=gene_name(existIndex);
%         sum(existIndex)
    
        log_data_zvalue=(log_data-repmat(mean(log_data,2),1,size(log_data,2)))./repmat(std(log_data')',1,size(log_data,2));
        log_data_doublets_zvalue=(log_data_doublets-repmat(mean(log_data_doublets,2),1,size(log_data_doublets,2)))./repmat(std(log_data_doublets')',1,size(log_data_doublets,2));
        log_data_zvalue(isnan(log_data_zvalue))=0;
        log_data_doublets_zvalue(isnan(log_data_doublets_zvalue))=0;
    
        folderName= strcat(n_g_input_path,'/data');
        foldername2= strcat(folderName,'/heatmapPIC_vs_AD')
        mkdir(foldername2)
        
        %%%%%%%%%%%%% 1) DEG by clustering 6 %%%%%%%%%%%%%
        clusterSize=max(clustering_color);
        [pvalue_total,fdr_total,logRatio_total,zvalue_total]=DEG_ranksum4cluster(clusterSize,log_data,clustering_color);
        disp('DEG ranksum Done!')
    %%
    %%%%%%%%%%%%% ---- Pipeline for clustering ---- %%%%%%%%%%%%%
    
        save([folderName,'/pvalue_fdr_logRatio_zvalue.mat'],'pvalue_total','fdr_total','logRatio_total','zvalue_total')
        save([folderName,'/use_parameters.mat'],'n_e_co','n_p_co','n_l_co') 

        %%%%%%%%%%%%% 2) PIC SVM classification %%%%%%%%%%%%%
        clusterSelect = n_clusters;
        load([folderName,'/pvalue_fdr_logRatio_zvalue.mat'])
        pCutoff = n_p_co; lrCutoff = n_l_co;
        seedNumber=1;randSize=10000;
        DEGnumber=5;
        [bestMatch,artificialDoubletsCombiUnique,tt,svm_score_t]=NicheSVM(pvalue_total,pCutoff,logRatio_total,lrCutoff,seedNumber,randSize,clustering_color,clusterSelect,clustering_name_unique,log_data_zvalue,log_data_doublets_zvalue,DEGnumber);
        save([folderName,'/SVM_bestMatch.mat'],'bestMatch','artificialDoubletsCombiUnique')
        drawHeatmap_DEG(folderName,clustering_color,pvalue_total,pCutoff,logRatio_total,lrCutoff,log_data_zvalue,gene_name,DEGnumber,clustering_name,clustering_name_unique)
    
        %%% Draw heatmap
        load([folderName,'/pvalue_fdr_logRatio_zvalue.mat'])
        clusterSelect = 1:13;
        clusterSelect = n_clusters;
        outputFile=[folderName,'/heatmap_PICSVM_top5DEG.pdf'];
        outputFile2=[folderName,'/heatmap_PICSVM_top5DEG.fig'];
        DEGnumber=5;
        max_t = max(svm_score_t');
        scoreIndex = find(max_t>-0.0001);
        drawHeatmap_PICSVM(outputFile,bestMatch(scoreIndex),artificialDoubletsCombiUnique,clusterSelect,clustering_name_unique,pvalue_total,pCutoff,logRatio_total,lrCutoff,log_data_doublets_zvalue(:,scoreIndex),gene_name,DEGnumber);
        drawHeatmap_PICSVM(outputFile2,bestMatch(scoreIndex),artificialDoubletsCombiUnique,clusterSelect,clustering_name_unique,pvalue_total,pCutoff,logRatio_total,lrCutoff,log_data_doublets_zvalue(:,scoreIndex),gene_name,DEGnumber);
        %outputFile = 'E:\Niche_SVM_GUI_Final\NicheSVM_GUI\heatmap_PICSVM_top5DEG.pdf';
        %drawHeatmap_PICSVM(outputFile,bestMatch,artificialDoubletsCombiUnique,clusterSelect,clustering_name_unique,pvalue_total,n_p_co,logRatio_total,n_l_co,log_data_doublets_zvalue,gene_name,DEGnumber);
    
        
        %%%%%%%%%%%%% 3) DEG between PICs and artificial doublets %%%%%%%%%%%%%
        clusterSelect = n_clusters;
        seedNumber=1;randSize=10000;
        minCell=5;
        [pvalue_total,fdr_total,logRatio_total]=DEG_PIC_vs_AD(bestMatch,seedNumber,randSize,clustering_color,clusterSelect,log_data_zvalue,clustering_name_unique,log_data_doublets_zvalue,minCell);
        save([folderName,'/pvalue_fdr_logRatio_PIC_vs_AD.mat'],'pvalue_total','fdr_total','logRatio_total')
    
        outputFolder=[folderName,'/heatmapPIC_vs_AD'];
        load([folderName,'/pvalue_fdr_logRatio_PIC_vs_AD.mat'])
        pvalue_totalPIC_AD=pvalue_total;
        logRatio_totalPIC_AD=logRatio_total;
        load([folderName,'/pvalue_fdr_logRatio_zvalue.mat'])
        pCutoff=n_p_co;lrCutoff = n_l_co;
        DEGnumber=5;DEGnumberPIC_AD=10;
        drawHeatmap_PIC_vs_AD(outputFolder,bestMatch,artificialDoubletsCombiUnique,clustering_color,clusterSelect,clustering_name_unique,log_data_zvalue,pvalue_total,pCutoff,logRatio_total,lrCutoff,log_data_doublets_zvalue,pvalue_totalPIC_AD,logRatio_totalPIC_AD,gene_name,DEGnumber,DEGnumberPIC_AD);
        drawHeatmap_DEG(folderName,clustering_color,pvalue_total,pCutoff,logRatio_total,lrCutoff,log_data_zvalue,gene_name,DEGnumber,clustering_name,clustering_name_unique)
        
        mk_deglist = strcat(folderName,'/deglist')
        mkdir(mk_deglist)
        DEGlists_PIC_vs_AD(mk_deglist,bestMatch,artificialDoubletsCombiUnique,pCutoff,lrCutoff,pvalue_totalPIC_AD,logRatio_totalPIC_AD,gene_name)
        msgbox("Analysis Completed!","Success");



%         mail = 'jahanzeb.ssu@gmail.com'; 
%         password = 'babajee123';
%         setpref('Internet','SMTP_Server','smtp.gmail.com');
%         setpref('Internet','E_mail',mail);
%         setpref('Internet','SMTP_Username',mail);
%         setpref('Internet','SMTP_Password',password);
%         props = java.lang.System.getProperties;
%         props.setProperty('mail.smtp.auth','true');
%         props.setProperty('mail.smtp.starttls.enable','true');
%         sendmail('bscs133301@gmail.com','R SVM compilation completed','Hello! This is a test from MATLAB!')
        
    catch me
        disp('Something bad happend!')
%         mail = 'jahanzeb.ssu@gmail.com'; 
%         password = 'babajee123';
%         setpref('Internet','SMTP_Server','smtp.gmail.com');
%         setpref('Internet','E_mail',mail);
%         setpref('Internet','SMTP_Username',mail);
%         setpref('Internet','SMTP_Password',password);
%         props = java.lang.System.getProperties;
%         props.setProperty('mail.smtp.auth','true');
%         props.setProperty('mail.smtp.starttls.enable','true');
%         sendmail('bscs133301@gmail.com','R SVM compilation FAIL','Hello! This is a test from MATLAB!')
%         msgbox("Something went wrong!","Error","error");
        
    end
end

