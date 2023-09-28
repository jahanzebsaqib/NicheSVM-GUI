function preprocess_2(inputArg1)
%PREPROCESS_2 Summary of this function goes here
%   Detailed explanation goes here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(strcat(inputArg1,'\process_data\data_output.mat'))

log_data          = log2(rna_matrix + 1);
log_data_doublets = log2(st_matrix  + 1);

clear rna_matrix st_matrix

clustering_name_unique = unique(clustering_name);
clustering_color = zeros(length(clustering_name),1);

for i=1:length(clustering_name_unique)
    temp_cluster = clustering_name_unique(i);
    temp_index = find(strcmp(temp_cluster,clustering_name));
    clustering_color(temp_index,1) = i;
end

clustering_color = clustering_color';

clearvars -except inputArg1 gene_name clustering_color clustering_name clustering_name_unique log_data log_data_doublets st_cordinates

save_name = strcat(inputArg1,'\process_data\ready_vector.mat');
save(save_name,'-v7.3')

close all

disp('Data Preprocessing Completed!')
msgbox("Data Preprocessing Completed!","Success");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

