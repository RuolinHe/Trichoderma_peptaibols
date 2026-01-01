function newset=BGC_AddidforNPRSs_HRL20231220(oldset,rawset,regions_flag,omains_flag,CDSs_flag,cand_clusters_flag)
%% Update
% use xxx_flag to control whether read xxx information (xxx=regions,
% omains, CDSs)
%   By Ruolin He in 2022/10/30
%   Updated by Ruolin He in 2023/12/20
%%
if nargin<3
    regions_flag=1;
    omains_flag=1;
    CDSs_flag=1;
    cand_clusters_flag=1;
end
% one set include: genomes, regions, CDSs, omains, interomains
if regions_flag
    [regions,oldset,regions_fields] = initialize_check(rawset,oldset,'regions');
end
if cand_clusters_flag % add cand_clusters
    [cand_clusters,oldset,cand_clusters_fields] = initialize_check(rawset,oldset,'cand_clusters');
end
if CDSs_flag % add CDS
    [CDSs,oldset,CDSs_fields] = initialize_check(rawset,oldset,'CDSs');
end
if omains_flag % add domain and interdomains
    [omains,oldset,omains_fields] = initialize_check(rawset,oldset,'omains');
end
if regions_flag
    % [regions,oldset,regions_fields] = initialize_check(rawset,oldset,'regions');
    if isfield(regions,'TTA')
        for j = 1:regions.num
            for i = 1:regions.TTA(j).num
                regions.TTA(j).CDSs_id=regions.TTA(j).CDSs_id+oldset.CDSs.num;
            end
        end
    end
    if oldset.regions.num ==0
        oldset.regions = regions;
    else
        oldset.regions.num=oldset.regions.num+regions.num;
        for j = 1:length(regions_fields)
            if ~strcmp(regions_fields{j},'num')
                oldset.regions.(regions_fields{j})=[oldset.regions.(regions_fields{j});regions.(regions_fields{j})];
            end
        end
    end
end
if cand_clusters_flag % add cand_clusters
    % [cand_clusters,oldset,cand_clusters_fields] = initialize_check(rawset,oldset,'cand_clusters');
    if oldset.cand_clusters.num == 0
        oldset.cand_clusters = cand_clusters;
    else
        oldset.cand_clusters.num=oldset.cand_clusters.num+cand_clusters.num;
        for j = 1:length(cand_clusters_fields)
            if strcmp(cand_clusters_fields{j},'region_ids')
                oldset.cand_clusters.(cand_clusters_fields{j})=[oldset.cand_clusters.(cand_clusters_fields{j});cand_clusters.(cand_clusters_fields{j})+oldset.regions.num-regions.num]; %这一行需要check
            elseif ~strcmp(cand_clusters_fields{j},'num')
                oldset.cand_clusters.(cand_clusters_fields{j})=[oldset.cand_clusters.(cand_clusters_fields{j});cand_clusters.(cand_clusters_fields{j})];
            end
        end
    end
end
if CDSs_flag % add CDS
    % [CDSs,oldset,CDSs_fields] = initialize_check(rawset,oldset,'CDSs');
    if oldset.CDSs.num == 0
        oldset.CDSs=CDSs;
    else
        oldset.CDSs.num=oldset.CDSs.num+CDSs.num;
        for j = 1:length(CDSs_fields)
            if strcmp(CDSs_fields{j},'region_ids')
                oldset.CDSs.(CDSs_fields{j})=[oldset.CDSs.(CDSs_fields{j});CDSs.(CDSs_fields{j})+oldset.regions.num-regions.num]; %这一行需要check
            elseif ~strcmp(CDSs_fields{j},'num')
                oldset.CDSs.(CDSs_fields{j})=[oldset.CDSs.(CDSs_fields{j});CDSs.(CDSs_fields{j})];
            end
        end
    end
end
if omains_flag % add domain and interdomains
    % [omains,oldset,omains_fields] = initialize_check(rawset,oldset,'omains');
    if isfield(omains,'domaintypelist')
        if ~isfield(oldset.omains,'domaintypelist')||isempty(oldset.omains.domaintypelist)
            oldset.omains.domaintypelist=omains.domaintypelist;
        else % check if omains are consistent
            for k=1:length(omains.domaintypelist)
                if ~strcmp(oldset.omains.domaintypelist{k},omains.domaintypelist{k})
                    perror('type list is not same')
                end
            end
        end
    end
    if oldset.omains.num==0
        oldset.omains=omains;
    else
        oldset.omains.num=oldset.omains.num+omains.num;
        for j = 1:length(omains_fields)
            if strcmp(omains_fields{j},'region_ids')
                oldset.omains.(omains_fields{j})=[oldset.omains.(omains_fields{j});omains.(omains_fields{j})+oldset.regions.num-regions.num];
            elseif strcmp(omains_fields{j},'CDS_ids')
                oldset.omains.(omains_fields{j})=[oldset.omains.(omains_fields{j});omains.(omains_fields{j})+oldset.CDSs.num-CDSs.num];
            elseif ~strcmp(omains_fields{j},'domaintypelist')&&~strcmp(omains_fields{j},'num')
                oldset.omains.(omains_fields{j})=[oldset.omains.(omains_fields{j});omains.(omains_fields{j})];
            end
        end
    end
end
if regions_flag&&isfield(oldset.regions,'product_after_merge')
    for i = 1:length(oldset.regions.product_after_merge)
        for j = 1:length(oldset.regions.product_after_merge(i).cand_cluster_ids)
            oldset.regions.product_after_merge(i).cand_cluster_ids{j}=oldset.regions.product_after_merge(i).cand_cluster_ids{j}+oldset.regions.num-regions.num;
        end
    end
end
newset=oldset;
end

function [omains,oldset,omains_fields] = initialize_check(rawset,oldset,struct_name)
omains=rawset.(struct_name);
omains_fields = fieldnames(omains);
if isempty(oldset)||~isfield(oldset,struct_name) % initiate if there is not a oldset yet
    oldset.(struct_name)=[];
    for j = 1:length(omains_fields)
        oldset.(struct_name).(omains_fields{j})=[];
    end
    oldset.(struct_name).num=0;
elseif oldset.(struct_name).num~=0 % compatibility
    old_omains_fields = fieldnames(oldset.(struct_name));
    if sum(ismember(old_omains_fields,omains_fields))~=length(omains_fields)
       warning(['The fields of old and new ',struct_name,' are inconsistent! Only keep the sharing filed'])
       omains_fields = omains_fields(ismember(omains_fields,old_omains_fields));
       oldset.(struct_name) = rmfield(oldset.(struct_name),omains_fields(~ismember(omains_fields,old_omains_fields)));
    end
end
end