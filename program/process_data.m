% update in 2022.10.27 by Ruolin He:  handle species (folders) which don't contain SM
% update in 2022.10.27 by Ruolin He:  generate each one by parallel
% generate data and database
% run it first
clear
clc
% warning off
antismashversion=7;

tic
sepstr='/';
toolboxpath='../tools_1025/';
if exist(toolboxpath,'dir')
    addpath(toolboxpath)
else
    error('Wrong path for MATLAB additional toolbox, please check\n')
end


datapaths='../data/antismash/galaxy/';
mat_outputpath='./output/matsummary/';
finalsetpath='./output/final.mat';

regions_flag=1; % only save regions
omains_flag=1;
CDSs_flag=1;
cand_clusters_flag=1;
printprogress=0;

if ~exist(mat_outputpath,'dir')
    mkdir(mat_outputpath);
    save_flage=1;
else
    save_flage=0;
end
save_flage=1;

% record the folder name, whether it has a zip folder, and how many gbk
% files they got
%%
non_SM_region = [];
non_SM_region_foldername = []; % record folders which don't contain any secondary metabolism BGCs detected by antiSMASH

finalset=[];
foldersinfolder=dir([datapaths,'GCA_*']);
%% 
for j=1:length(foldersinfolder) % check each antismash folder
    fname=foldersinfolder(j).name;
    fprintf('%d/%d\t%s\n',j,length(foldersinfolder),fname)
    folderpath=[datapaths,fname,sepstr];
    if length(dir(folderpath))>2 % pass empty folder
        % if there is a zip file, read the zip file, output a .mat  file (if not yet exist)
        matmat_outputpath=[mat_outputpath,num2str(j),'_',fname,'.mat'];
        if ~exist(matmat_outputpath,'file') % not yet has the output
            %%%%%%%%%%%%%%%%%%%% Read GBKs in this folder for information
            rawset=GBK_Read_Antismash_HRL(datapaths,fname,sepstr,printprogress,antismashversion,regions_flag,omains_flag,CDSs_flag,cand_clusters_flag);
            % save the raw set into mat
            save(matmat_outputpath,'rawset');
        else % already have the raw set, just load it.
            load(matmat_outputpath,'rawset');
        end
        finalset=BGC_AddidforNPRSs_HRL20231220(finalset,rawset,regions_flag,omains_flag,CDSs_flag,cand_clusters_flag);
        if rawset.regions.num == 0 % record species (folders) which don't contain SM
            non_SM_region_foldername=[non_SM_region_foldername;{[datapaths,fname]}];
        end
    else
        error('%s is empty!',fname)
    end
end
fname='Trichoderma_hypoxylon';
matmat_outputpath=[mat_outputpath,num2str(j+1),'_',fname,'.mat'];
if ~exist(matmat_outputpath,'file') % not yet has the output
    %%%%%%%%%%%%%%%%%%%% Read GBKs in this folder for information
    rawset=GBK_Read_Antismash_HRL('../data/antismash/',fname,sepstr,printprogress,antismashversion,regions_flag,omains_flag,CDSs_flag,cand_clusters_flag);
    % save the raw set into mat
    save(matmat_outputpath,'rawset');
else % already have the raw set, just load it.
    load(matmat_outputpath,'rawset');
end
finalset=BGC_AddidforNPRSs_HRL20231220(finalset,rawset,regions_flag,omains_flag,CDSs_flag,cand_clusters_flag);
%% 
if ~isempty(non_SM_region_foldername)
    non_SM_region.foldername=non_SM_region_foldername;
end
%%
if save_flage % only save once (first time)
    fprintf('saving data...');
    my_regions=finalset.regions;
    save('./output/my_regions','my_regions');
    my_cand_clusters=finalset.cand_clusters;
    save('./output/my_cand_clusters','my_cand_clusters');
    my_omains=finalset.omains;
    save('./output/my_omains','my_omains');
    my_CDSs=finalset.CDSs;
    save('./output/my_CDSs','my_CDSs','-v7.3');
    save('./output/non_SM_region','non_SM_region');
    fprintf('done')
end
toc
