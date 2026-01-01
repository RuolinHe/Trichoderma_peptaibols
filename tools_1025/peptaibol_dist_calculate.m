function best_peptaibol = peptaibol_dist_calculate(best_peptaibol,ref_index,motif_sepcial_align)
%peptaibol_dist_calculate Summary of this function goes here
% input:
% ref_index=find(ismember(best_peptaibol.peptaibol_name,'LogPS1'));
% 20 module peptaibol as reference
% motif_sepcial_align: n*2 cell. The first is orginal motif, the second is
% aligned motif (maunally). For some motifs which don't have same length
% output
%   msa_list
%     A3_A6_dist_matrix: p-distance
%     A3_A6_align_module
%     A3_A6_align_module_min_dist
%     A3_A6_align_module_dist_refer
%     A3_A6_min_module
%     A3_A6_min_module_dist
%     A3_A6_same_order: for each peptaibol, if they have same module order
%     A3_A6_merged_group: merge peptaibol with same module order into group
%     A3_A6_dist_matrix_within_PS: dist within peptaibol synthatse
%     A4_A5_dist_matrix
%     A4_A5_align_module
%     A4_A5_align_module_min_dist
%     A4_A5_align_module_dist_refer
%     A4_A5_min_module
%     A4_A5_min_module_dist
%     A4_A5_same_order
%     A4_A5_merged_group
%     A4_A5_dist_matrix_within_PS
%% 
ref_id_index=find(best_peptaibol.peptaibol_list==ref_index);
loc_pwd=pwd;
if exist('tmp_24601','dir')==7
    rmdir('tmp_24601', 's')
end
mkdir tmp_24601 % bulid a temporary folder
best_peptaibol.msa_list=cell(1,size(best_peptaibol.seq_list,2));
for k = 1:size(best_peptaibol.seq_list,2)
% for k = size(best_peptaibol.seq_list,2)
    if mod(k,2)~=0%motif
        best_peptaibol.msa_list{k}=[];
        for i = 1:size(best_peptaibol.seq_list,1)
            loc_index=ismember(motif_sepcial_align(:,1),best_peptaibol.seq_list{i,k});
            if any(loc_index)
                best_peptaibol.msa_list{k}=[best_peptaibol.msa_list{k};motif_sepcial_align{loc_index,2}];
            else
                best_peptaibol.msa_list{k}=[best_peptaibol.msa_list{k};best_peptaibol.seq_list{i,k}];
            end
        end
    else%intermotif
        test_header = cellstr(num2str([1:size(best_peptaibol.seq_list,1)]'));
        fastawrite([pwd,'/tmp_24601/',num2str(k),'.fasta'],test_header,best_peptaibol.seq_list(:,k));
%         command = ['clustalo -i ',loc_pwd,'/tmp_24601/',num2str(k),'.fasta -o ',loc_pwd,'/tmp_24601/',num2str(k),'_1.fasta --outfmt fasta --output-order input-order'];
        command = ['muscle -super5 ',loc_pwd,'/tmp_24601/',num2str(k),'.fasta -output ',loc_pwd,'/tmp_24601/',num2str(k),'_1.fasta'];% too much sequences, use muscle rather than clustalo
        system(command);
        best_peptaibol.msa_list{k}=[];
        raw = fastaread([pwd,'/tmp_24601/',num2str(k),'_1.fasta']);
        msa_order=[];
        loc_msa=[];
        for ii = 1:length(raw)
            msa_order=[msa_order;str2double(raw(ii).Header)];
            loc_msa=[loc_msa;raw(ii).Sequence];
        end
        [~,I]=sort(msa_order,'ascend');
        loc_msa=loc_msa(I,:);
        best_peptaibol.msa_list{k}=loc_msa;
    end
end
rmdir('tmp_24601', 's')
%%
loc_seq=[];
for i = 1:length(best_peptaibol.msa_list)
    loc_seq=[loc_seq,best_peptaibol.msa_list{i}];
end
% best_peptaibol.A3_A6_dist_matrix=seqpdist(loc_seq,'ScoringMatrix','BLOSUM62','Method','alignment-score','SquareForm',true);
best_peptaibol.A3_A6_dist_matrix=seqpdist(loc_seq,'ScoringMatrix','BLOSUM62','Method','p-distance','SquareForm',true);
%% 
% take LogPS1 as reference because it has 20 modules and the most number of
% products
% ref_index=find(best_peptaibol.peptaibol_list==find(ismember(best_peptaibol.peptaibol_name,'LogPS1')));
best_peptaibol.A3_A6_align_module=nan(length(best_peptaibol.module_list),length(ref_id_index));
best_peptaibol.A3_A6_align_module_min_dist=nan(length(best_peptaibol.module_list),length(ref_id_index));
best_peptaibol.A3_A6_align_module_dist_refer=[];
for i = 1:best_peptaibol.num
    loc_index=find(best_peptaibol.peptaibol_list==i);
    best_peptaibol.A3_A6_align_module_dist_refer=[best_peptaibol.A3_A6_align_module_dist_refer;best_peptaibol.A3_A6_dist_matrix(loc_index,ref_id_index)];
    for j = 1:length(loc_index)
        loc_dist=best_peptaibol.A3_A6_dist_matrix(loc_index(j),ref_id_index);
        [min_dist,I]=min(loc_dist);
        best_peptaibol.A3_A6_align_module(loc_index(j),I)=j;
        best_peptaibol.A3_A6_align_module_min_dist(loc_index(j),I)=min_dist;
    end
end
%% 
best_peptaibol.A3_A6_align_module_matrix=cell(best_peptaibol.num,length(ref_id_index));
for i = 1:best_peptaibol.num
    loc_align_module=best_peptaibol.A3_A6_align_module(best_peptaibol.peptaibol_list==i,:);
    for j = 1:size(best_peptaibol.A3_A6_align_module,2)
        loc_module=loc_align_module(~isnan(loc_align_module(:,j)),j);
        if ~isempty(loc_module)
            if length(loc_module)==1
                best_peptaibol.A3_A6_align_module_matrix(i,j)=cellstr(num2str(loc_module));
            else
                best_peptaibol.A3_A6_align_module_matrix(i,j)=join(cellfun(@num2str,num2cell(loc_module),'UniformOutput',false),'/');
            end
        end
    end
end
%% 
best_peptaibol.A3_A6_min_module=nan(length(best_peptaibol.peptaibol_list),best_peptaibol.num);
best_peptaibol.A3_A6_min_module_dist=nan(length(best_peptaibol.peptaibol_list),best_peptaibol.num);
for i = 1:length(best_peptaibol.peptaibol_list)
    for j = 1:best_peptaibol.num
        [min_dist,I]=min(best_peptaibol.A3_A6_dist_matrix(i,best_peptaibol.peptaibol_list==j));
        best_peptaibol.A3_A6_min_module(i,j)=I;
        best_peptaibol.A3_A6_min_module_dist(i,j)=min_dist;
    end
end
%% 
best_peptaibol.A3_A6_same_order=cell(best_peptaibol.num,1);
for i = 1:best_peptaibol.num
    for j = 1:best_peptaibol.num
        if all(ismember(best_peptaibol.A3_A6_min_module(best_peptaibol.peptaibol_list==i,i),best_peptaibol.A3_A6_min_module(best_peptaibol.peptaibol_list==i,j)))
            best_peptaibol.A3_A6_same_order{i}=[best_peptaibol.A3_A6_same_order{i},j];
        end
    end
end
%% 
merge_flags = true(1, numel(best_peptaibol.A3_A6_same_order));
for i = 1:numel(best_peptaibol.A3_A6_same_order)
    if merge_flags(i)
        for j = i+1:numel(best_peptaibol.A3_A6_same_order)
            % 检查数组之间是否有交集
            if ~isempty(intersect(best_peptaibol.A3_A6_same_order{i}, best_peptaibol.A3_A6_same_order{j}))
                % 合并数组并更新标志向量
                best_peptaibol.A3_A6_same_order{i} = union(best_peptaibol.A3_A6_same_order{i}, best_peptaibol.A3_A6_same_order{j});
                merge_flags(j) = false;
            end
        end
    end
end
best_peptaibol.A3_A6_merged_group = best_peptaibol.A3_A6_same_order(merge_flags);
%% 
best_peptaibol.A3_A6_dist_matrix_within_PS=cell(best_peptaibol.num,1);
for i = 1:best_peptaibol.num
    best_peptaibol.A3_A6_dist_matrix_within_PS{i}=best_peptaibol.A3_A6_dist_matrix(best_peptaibol.peptaibol_list==i,best_peptaibol.peptaibol_list==i);
end
%%
loc_seq=[];
for i = 3:5
    loc_seq=[loc_seq,best_peptaibol.msa_list{i}];
end
% best_peptaibol.A3_A6_dist_matrix=seqpdist(loc_seq,'ScoringMatrix','BLOSUM62','Method','alignment-score','SquareForm',true);
best_peptaibol.A4_A5_dist_matrix=seqpdist(loc_seq,'ScoringMatrix','BLOSUM62','Method','p-distance','SquareForm',true);
%% 
% take LogPS1 as reference because it has 20 modules and the most number of
% products
% ref_index=find(best_peptaibol.peptaibol_list==find(ismember(best_peptaibol.peptaibol_name,'LogPS1')));
best_peptaibol.A4_A5_align_module=nan(length(best_peptaibol.module_list),length(ref_id_index));
best_peptaibol.A4_A5_align_module_min_dist=nan(length(best_peptaibol.module_list),length(ref_id_index));
best_peptaibol.A4_A5_align_module_dist_refer=[];
for i = 1:best_peptaibol.num
    loc_index=find(best_peptaibol.peptaibol_list==i);
    best_peptaibol.A4_A5_align_module_dist_refer=[best_peptaibol.A4_A5_align_module_dist_refer;best_peptaibol.A4_A5_dist_matrix(loc_index,ref_id_index)];
    for j = 1:length(loc_index)
        loc_dist=best_peptaibol.A4_A5_dist_matrix(loc_index(j),ref_id_index);
        [min_dist,I]=min(loc_dist);
        best_peptaibol.A4_A5_align_module(loc_index(j),I)=j;
        best_peptaibol.A4_A5_align_module_min_dist(loc_index(j),I)=min_dist;
    end
end
%% 
best_peptaibol.A4_A5_align_module_matrix=cell(best_peptaibol.num,length(ref_id_index));
for i = 1:best_peptaibol.num
    loc_align_module=best_peptaibol.A4_A5_align_module(best_peptaibol.peptaibol_list==i,:);
    for j = 1:size(best_peptaibol.A4_A5_align_module,2)
        loc_module=loc_align_module(~isnan(loc_align_module(:,j)),j);
        if ~isempty(loc_module)
            if length(loc_module)==1
                best_peptaibol.A4_A5_align_module_matrix(i,j)=cellstr(num2str(loc_module));
            else
                best_peptaibol.A4_A5_align_module_matrix(i,j)=join(cellfun(@num2str,num2cell(loc_module),'UniformOutput',false),'/');
            end
        end
    end
end
%%
best_peptaibol.A4_A5_min_module=nan(length(best_peptaibol.peptaibol_list),best_peptaibol.num);
best_peptaibol.A4_A5_min_module_dist=nan(length(best_peptaibol.peptaibol_list),best_peptaibol.num);
for i = 1:length(best_peptaibol.peptaibol_list)
    for j = 1:best_peptaibol.num
        [min_dist,I]=min(best_peptaibol.A4_A5_dist_matrix(i,best_peptaibol.peptaibol_list==j));
        best_peptaibol.A4_A5_min_module(i,j)=I;
        best_peptaibol.A4_A5_min_module_dist(i,j)=min_dist;
    end
end
%% 
best_peptaibol.A4_A5_same_order=cell(best_peptaibol.num,1);
for i = 1:best_peptaibol.num
    for j = 1:best_peptaibol.num
        if all(ismember(best_peptaibol.A4_A5_min_module(best_peptaibol.peptaibol_list==i,i),best_peptaibol.A4_A5_min_module(best_peptaibol.peptaibol_list==i,j)))
            best_peptaibol.A4_A5_same_order{i}=[best_peptaibol.A4_A5_same_order{i},j];
        end
    end
end
%% 
merge_flags = true(1, numel(best_peptaibol.A4_A5_same_order));
for i = 1:numel(best_peptaibol.A4_A5_same_order)
    if merge_flags(i)
        for j = i+1:numel(best_peptaibol.A4_A5_same_order)
            if ~isempty(intersect(best_peptaibol.A4_A5_same_order{i}, best_peptaibol.A4_A5_same_order{j}))
                best_peptaibol.A4_A5_same_order{i} = union(best_peptaibol.A4_A5_same_order{i}, best_peptaibol.A4_A5_same_order{j});
                merge_flags(j) = false;
            end
        end
    end
end
best_peptaibol.A4_A5_merged_group = best_peptaibol.A4_A5_same_order(merge_flags);
%% 
best_peptaibol.A4_A5_dist_matrix_within_PS=cell(best_peptaibol.num,1);
for i = 1:best_peptaibol.num
    best_peptaibol.A4_A5_dist_matrix_within_PS{i}=best_peptaibol.A4_A5_dist_matrix(best_peptaibol.peptaibol_list==i,best_peptaibol.peptaibol_list==i);
end
end

