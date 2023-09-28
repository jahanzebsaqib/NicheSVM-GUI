function preprocess_1(inputArg1,inputArg2,inputArg3)
%PREPROCESS_1 Summary of this function goes here
%   Detailed explanation goes here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath(genpath('C:\Users\ja_ha\OneDrive\Documents\MATLAB\scRNA_ST_intergration_PDAC\PDAC-B\Chapter_1'));

disp('Preprocessing Started!')

rna_table = readtable(inputArg1,'Format','auto','ReadVariableNames',false,'ReadRowNames',false);
clustering_name = string(table2cell(rna_table(1,2:end)));

disp('scRNA file loaded!')


clear rna_table

rna_table = readtable(inputArg1,'ReadVariableNames',false,'ReadRowNames',false);
rna_gene_names = string(table2cell(rna_table(:,1)));
raw_rna_data = table2array(rna_table(:,2:end));

clear rna_table

st_table = readtable(inputArg2,'Format','auto','FileType','text','ReadVariableNames',false,'ReadRowNames',false);
st_cordinates = string(table2cell(st_table(1,2:end)));
disp('ST file loaded!')

clear st_table

st_table = readtable(inputArg2,'FileType','text','ReadVariableNames',false,'ReadRowNames',false);
st_gene_names = string(table2cell(st_table(:,1)));
raw_st_data = table2array(st_table(:,2:end));

clear st_table

gene_name = intersect(rna_gene_names,st_gene_names, 'stable');

%%
rna_matrix = zeros(length(gene_name(:,1)),length(raw_rna_data(1,:)));
st_matrix  = zeros(length(gene_name(:,1)),length(raw_st_data(1,:)));

for i=1:length(gene_name(:,1))
    temp_gene = gene_name(i);
    rna_index = find(strcmp(temp_gene,rna_gene_names));
    st_index  = find(strcmp(temp_gene,st_gene_names));
    rna_matrix(i,:) = ceil(mean(raw_rna_data(rna_index,:),1));
    st_matrix(i,:)  = ceil(mean(raw_st_data(st_index,:),1));
end

clearvars -except rna_matrix st_matrix st_cordinates clustering_name gene_name inputArg3

mk_name = strcat(inputArg3,'\process_data');
mkdir(mk_name)

save_name = strcat(inputArg3,'\process_data\data_output.mat');
save(save_name,'-v7.3')

close all
clc

disp('Initial Data Preprocessing Completed!')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

