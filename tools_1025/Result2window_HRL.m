function [aa_window,nt_window,my_start,my_end,motifid_mat,index_mat] = Result2window_HRL(result,NT_seq,Q_start,Q_end,window,slide)
%Result2window_HRL 把Find_NRPS_motif_module_pfam_HRL的结果转化成滑窗结果
%input
%   result:Find_NRPS_motif_module_pfam_HRL的结果
%   NT_seq可以为空则nt_window为空
%   Q_start：从哪个开始，第一个数字是第几个，第二个数字是什么domain.C=1,A=2,T=3;
%   Q_end：到哪个结束（包含该module) module=第一个motif到最后一个motif和最后一个后面的intermotif
%   window:滑窗长度
%   slide：滑窗步长
%output
%   aa_window,nt_window:滑窗结果
%   my_start：第一个起始的氨基酸的索引
%   my_end：最后一个氨基酸的索引
%   motifid_mat：2列，只含motif的位置
%   index_mat：motifid_mat的氨基酸起始索引
start_index=find(ismember(result.motifid_mat,[Q_start(2),1],'rows'));
end_index=find(ismember(result.motifid_mat,[Q_end(2),1],'rows'));
loc_end_index = end_index(Q_end(1));
end_domain=result.motifid_mat(loc_end_index,1);
loc_domain=end_domain;
while loc_end_index+1<=length(result.motifid_mat)&&loc_domain == end_domain
    loc_end_index=loc_end_index+1;
    loc_domain=result.motifid_mat(loc_end_index,1);
end
aa_seq=[];
for i = start_index(Q_start(1)):loc_end_index
    aa_seq=[aa_seq,result.seq_list{i}];
end
aa_window = Seq2window_HRL(aa_seq,window,slide);
my_start=result.index_mat(start_index(Q_start(1)),1);
my_end=result.index_mat(loc_end_index,2);
if ~isempty(NT_seq)
    nt_seq=NT_seq(3*my_start-2:3*my_end);
    nt_window = Seq2window_HRL(nt_seq,window*3,slide*3);
else
    nt_window=[];
end
motifid_mat = result.motifid_mat(start_index(Q_start(1)):loc_end_index,:);
index_mat = result.index_mat(start_index(Q_start(1)):loc_end_index,:);
my_rem = rem(motifid_mat(:,2),1);
motifid_mat(my_rem~=0,:)=[];
index_mat(my_rem~=0,:)=[];
end

