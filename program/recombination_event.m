toolboxpath='../tools_1025/';
if exist(toolboxpath,'dir')
    addpath(toolboxpath)
else
    error('Wrong path for MATLAB additional toolbox, please check\n')
end
if ~exist('best_peptaibol','var')
    load("Known_peptaibol_struct.mat")
end
%%
% DNA, window=150nt, step=30nt
CAT_color=[0,112,224;249,126,0;226,226,0]/255;
aa_window=50;
aa_slide=10;
%% 32 vs 9
% 32, M1-2, M17-18 (shorter is query)
% 9, M1-3, M18-20 (longer is reference)
event_pair_list={[32,1,2,3;9,1,3,3];[32,17,18,3;9,18,20,3];[18,4,7,1;23,4,7,1];[24,16,18,3;25,16,18,3];[19,4,5,3;27,4,6,3];[16,5,6,3;27,5,7,3];[6,8,10,3;16,8,10,3];[36,16,17,3;6,16,18,3];[5,16,17,3;6,16,18,3];[13,4,6,3;6,4,6,3];[12,4,6,3;6,4,6,3];[1,5,7,3;1,16,18,3];[9,1,3,3;9,13,15,3];[13,4,6,3;13,11,13,3]};
%% 
for event_index = 1:length(event_pair_list)
% for event_index = 1:2
    event_pair=event_pair_list{event_index};
    % event_pair=[32,1,2;9,1,3];%[query_index,query_start_module_index,query_end_module_index,end_domain_id;reference_index,reference_start_module_index,reference_end_module_index,end_domain_id];
    if event_pair(1,4)==3 % Before NRPS, there is an additional T domain in PKS
        [Q_aa_window,Q_nt_window,Q_start,Q_end,Q_motifid_mat,Q_index_mat] = Result2window_HRL(best_long_peptaibol.result_list{event_pair(1,1)},best_long_peptaibol.NT_seqs_list{event_pair(1,1)},[event_pair(1,2),1],[event_pair(1,3)+1,event_pair(1,4)],aa_window,aa_slide);
    else
        [Q_aa_window,Q_nt_window,Q_start,Q_end,Q_motifid_mat,Q_index_mat] = Result2window_HRL(best_long_peptaibol.result_list{event_pair(1,1)},best_long_peptaibol.NT_seqs_list{event_pair(1,1)},[event_pair(1,2),1],[event_pair(1,3),event_pair(1,4)],aa_window,aa_slide);
    end
    if event_pair(2,4)==3
        [R_aa_window,R_nt_window,R_start,R_end,R_motifid_mat,R_index_mat] = Result2window_HRL(best_long_peptaibol.result_list{event_pair(2,1)},best_long_peptaibol.NT_seqs_list{event_pair(2,1)},[event_pair(2,2),1],[event_pair(2,3)+1,event_pair(2,4)],aa_window,aa_slide);
    else
        [R_aa_window,R_nt_window,R_start,R_end,R_motifid_mat,R_index_mat] = Result2window_HRL(best_long_peptaibol.result_list{event_pair(2,1)},best_long_peptaibol.NT_seqs_list{event_pair(2,1)},[event_pair(2,2),1],[event_pair(2,3),event_pair(2,4)],aa_window,aa_slide);
    end
     
    score_matrix=zeros(length(Q_aa_window),length(R_aa_window));
    for j = 1:length(Q_aa_window)
        for i = 1:length(R_aa_window)
            [~,aln]=nwalign(Q_aa_window{j},R_aa_window{i},'ScoringMatrix','BLOSUM62');
            loc_NT_MSA=[seqinsertgaps(Q_nt_window{j},aln(1,:));seqinsertgaps(R_nt_window{i},aln(3,:))];
            loc_score_matrix=seqpdist(loc_NT_MSA,'ScoringMatrix','NUC44','Method','p-distance','SquareForm',true,'Alphabet','NT');
            score_matrix(j,i)=loc_score_matrix(1,2);
        end
    end
    score_matrix=1-score_matrix;
    search_range=30;
    [max_score,I]=max(score_matrix);
    [max_score1,I1]=max(score_matrix,[],2);
    Q_strain=best_long_peptaibol.Strain_used{event_pair(1,1)};
    R_strain=best_long_peptaibol.Strain_used{event_pair(2,1)};
    Q_strain=strrep(Q_strain,'Trichoderma','T.');
    R_strain=strrep(R_strain,'Trichoderma','T.');
    Q_str=['module ',num2str(event_pair(1,2)),'-',num2str(event_pair(1,3)),' from ',Q_strain];
    R_str=['module ',num2str(event_pair(2,2)),'-',num2str(event_pair(2,3)),' from ',R_strain];
    F_window_index_R=(R_index_mat*3)/1000;
    R_position=(((1:length(R_nt_window))-1)*30+R_start*3-2)/1000;
    F_window_index_Q=(Q_index_mat*3)/1000;
    Q_position=(((1:length(Q_nt_window))-1)*30+Q_start*3-2)/1000;
     
    figure('Units','normalized','outerposition',[0 0 1 1]);
    if event_pair(2,3)-event_pair(2,2) + 1>3
        subplot(2,3,1:3)
    else
        subplot(2,3,1:(event_pair(2,3)-event_pair(2,2)+1))
    end
    yyaxis right
    plot(R_position,((I-1)*30+Q_start*3-2)/1000,'g')
    loc_ylabel=['Position of maximum in ',Q_str,' (kb)'];
    loc_ylabel=split(loc_ylabel,' from');
    loc_ylabel{2}=['from ',loc_ylabel{2}];
    ylabel(loc_ylabel)
    yyaxis left
    hold on
    for i=1:size(F_window_index_R,1)
        patch([F_window_index_R(i,1),F_window_index_R(i,2),F_window_index_R(i,2),F_window_index_R(i,1)],[0 0 1 1],CAT_color(R_motifid_mat(i,1),:),'EdgeColor','none'); 
    end
    plot(R_position,max_score,'k')
    xlabel(['Position in ',R_str,' (kb)'])
    ylabel('Maximum sequence identity')
    ylim([0,1])
    hold off
    set(gca, 'Fontname', 'Arial');
    set(gca, 'Fontsize', 21);
    if event_pair(1,3)-event_pair(1,2)+1 >3
        subplot(2,3,4:6)
    else
        subplot(2,3,(1:(event_pair(1,3)-event_pair(1,2)+1))+3)
    end
    yyaxis right
    plot(Q_position,((I1-1)*30+R_start*3-2)/1000,'g')
    loc_ylabel=['Position of maximum in ',R_str,' (kb)'];
    loc_ylabel=split(loc_ylabel,' from');
    loc_ylabel{2}=['from ',loc_ylabel{2}];
    ylabel(loc_ylabel)
    yyaxis left
    hold on
    for i=1:size(F_window_index_Q,1)
        patch([F_window_index_Q(i,1),F_window_index_Q(i,2),F_window_index_Q(i,2),F_window_index_Q(i,1)],[0 0 1 1],CAT_color(Q_motifid_mat(i,1),:),'EdgeColor','none'); 
    end
    plot(Q_position,max_score1,'k')
    xlabel(['Position in ',Q_str,' (kb)'])
    ylabel('Maximum sequence identity')
    ylim([0,1])
    hold off
    set(gca, 'Fontname', 'Arial');
    set(gca, 'Fontsize', 21);
    saveas(gcf,['./output/figure/recombination_event/Event_',num2str(event_index),'.svg'])
end