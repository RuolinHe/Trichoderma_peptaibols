%% generate peptaibol struct
% from line 1-1454, we generate best_peptaibol, best_long_peptaibol,
% good_peptaibol, long_peptaibol and known_peptaibol_struct
% and these variables are saved in "known_peptaibol_struct.mat"
%
% form line 1455-1479, prepare for input of iTOL for A domain in Moudule 1
%
% form line 1480-1727, prepare for input of iTOL for Trichoderma speices
%
% from line 1728-1745, prepare for input of MLGO
%
% from line 1746-1812, UMAP based on sequence distance of A domains recognizing 
% Ala, Aib, Vxx, or Lxx
%%
sepstr='/';
toolboxpath='../tools_1025/';
if exist(toolboxpath,'dir')
    addpath(toolboxpath)
else
    error('Wrong path for MATLAB additional toolbox, please check\n')
end
%% 
if ~exist('my_omains','var')
    tic
    fprintf('Loading data...\t')
    load('./output/my_omains.mat');
    load('./output/my_regions.mat')
    load('./output/non_SM_region.mat')
    fprintf('Finsh\n')
    toc
end
my_regions.uni_sorted_product_type_str=cell(my_regions.num,1);
for i = 1:my_regions.num
    my_regions.uni_sorted_product_type_str(i)=join(my_regions.product_single_type{i},'+');
end
if ~exist('./output/figure','dir')
    mkdir('./output/figure');
end
A_domain=2;
%% 
fprintf('Find A domain motif...\t')
if ~exist('./output/result_T.mat','file')
    A_domain_index=find(my_omains.typeid_mat(:,1)==A_domain&my_omains.isdomain==1);
    result=cell(length(A_domain_index),1);
    for i = 1:length(A_domain_index)
        seq=[];
        seq.Header=['omain_ids|',num2str(A_domain_index(i))];
        seq.Sequence=[my_omains.seq_ntaa{A_domain_index(i)-1,2},my_omains.seq_ntaa{A_domain_index(i),2},my_omains.seq_ntaa{A_domain_index(i)+1,2}];
        if ~isempty(seq.Sequence)
            result{i} = Find_NRPS_motif_module_pfam_HRL(seq,[],[],{'Aalpha','G','Talpha'},0,[],[]);
        end
    end
    fprintf('Finsh\n')
    fprintf('Saving data...\t')
    save('./output/result_T','result','A_domain_index','-v7.3')
else
    load('./output/result_T.mat')
end
fprintf('Finsh\n')
%% 
tic
fprintf('Find Aib...\t')
msa_outputpath='./output/msasummary_T/';
if ~exist(msa_outputpath,'dir')
    mkdir(msa_outputpath);
    reference=fastaread('Aib_fix_A4_A5.fasta');
    for i = 1:length(result)
        result{i}.A_domain_index=A_domain_index(i);
        result{i}.good=0;
        if isfield(result{i},'domain_list')&&any(result{i}.domain_list==A_domain)
           [~,loc_index]=ismember([A_domain,5],result{i}.motifid_mat,'rows');%A4
           if ~isempty(result{i}.seq_list{loc_index})&&~isempty(result{i}.seq_list{loc_index+1})&&~isempty(result{i}.seq_list{loc_index+2}) % A4-A5 aren't empty
                result{i}.good=1;
                if ~exist([msa_outputpath,num2str(i),'_seq_MSA.fasta'],'file')
                    query=[];
                    query.Sequence=[result{i}.seq_list{loc_index},result{i}.seq_list{loc_index+1},result{i}.seq_list{loc_index+2}];
                    query.Header='query';
                    query=[query;reference];
                    fastawrite([msa_outputpath,num2str(i),'_seq.fasta'],query);
                    command = ['clustalo -i ',msa_outputpath,num2str(i),'_seq.fasta -o ',msa_outputpath,num2str(i),'_seq_MSA.fasta --outfmt fasta --output-order input-order --threads=64'];
                    system(command);
                end
                loc_seq_MSA=fastaread([msa_outputpath,num2str(i),'_seq_MSA.fasta']);
                seq_MSA=[];
                for j = 1:length(loc_seq_MSA)
                    seq_MSA=[seq_MSA;loc_seq_MSA(j).Sequence];
                end
                distm= seqpdist(seq_MSA,'ScoringMatrix','BLOSUM62','Method','alignment-score','SquareForm',true);
                result{i}.MinDist=min(distm(1,2:end));
                result{i}.MeanDist=mean(distm(1,2:end));
           end
        end
    end
    toc
    fprintf('Finsh\n')
    fprintf('Saving data...\t')
    save('./output/result_T','result','A_domain_index','-v7.3')
end
fprintf('Finsh\n')
%% 
Dist_matrix=NaN(length(result),2);
for i = 1:length(result)
    if result{i}.good==1
        Dist_matrix(i,1)=result{i}.MinDist;
        Dist_matrix(i,2)=result{i}.MeanDist;
    end
end
%% 
sum(Dist_matrix(:,1)<0.075) % number of Aib-A in Trichoderma
%% 
fprintf('generating folders...\t')
for r=1:my_regions.num
    flist=strsplit(my_regions.foldername{r},'/');
    my_regions.shortfolder{r,1}=flist{end};
end
all_shortfolder=my_regions.shortfolder;
if ~isempty(non_SM_region) % handle species (folders) which don't contain SM
    for r = 1:length(non_SM_region.foldername)
        flist=strsplit(non_SM_region.foldername{r},'/');
        non_SM_region.shortfolder{r,1}=flist{end};
    end
    all_shortfolder=[all_shortfolder;non_SM_region.shortfolder];
end

folders=[];
[folders.foldername,~,folderid_all]=unique(all_shortfolder);

tabfreq=tabulate(folderid_all);
folders.num=length(folders.foldername);
folders.region_nums=zeros(folders.num,1);
folders.region_nums(tabfreq(:,1))=tabfreq(:,2);
if ~isempty(non_SM_region)
    folders.region_nums(ismember(folders.foldername,non_SM_region.shortfolder))=zeros(sum(ismember(folders.foldername,non_SM_region.shortfolder)),1);
end
folderid_contain_SM = folderid_all;
folderid_contain_SM(folders.region_nums(folderid_contain_SM)==0)=[];
my_regions.folderid=folderid_contain_SM;
%% 
Assembly=readcell('../data/genome/fasta/Trichoderma_fasta.xlsx','Sheet','Genome Assembly Data Report','Range','A:A');
Species=readcell('../data/genome/fasta/Trichoderma_fasta.xlsx','Sheet','Genome Assembly Data Report','Range','O:O');
Known_peptaibol=readcell('../data/genome/fasta/Trichoderma_fasta.xlsx','Sheet','Genome Assembly Data Report','Range','P:P');
OrganismName=readcell('../data/genome/fasta/Trichoderma_fasta.xlsx','Sheet','Genome Assembly Data Report','Range','C:C');
folders.Species=cell(folders.num,1);
folders.Short_species=cell(folders.num,1);
folders.Known_peptaibol=ones(folders.num,1);
for i = 1:folders.num
    if strcmp(folders.foldername{i},'Trichoderma_hypoxylon')
        folders.Species{i}='T. hypoxylon';
        folders.Short_species{i}='Trichoderma hypoxylon';
    else
        folders.Species{i}=Species{ismember(Assembly,folders.foldername{i})};
        loc_name=split(OrganismName{ismember(Assembly,folders.foldername{i})});
        folders.Short_species(i)=join(loc_name(1:2));
        folders.Known_peptaibol(i)=Known_peptaibol{ismember(Assembly,folders.foldername{i})};
    end
end
%% 
figure('Units','normalized','outerposition',[0 0 1 1]);
histogram(Dist_matrix(:,1))
xlabel('Min distance to Aib')
ylabel('Count')
title(['A domain min distance distribution in the ',num2str(folders.num),' Trichoderma genomes (n=',num2str(sum(~isnan(Dist_matrix(:,1)))),')'])
saveas(gcf,'./output/figure/A_min_dist_distribution_Trichoderma.svg')
%% 
figure('Units','normalized','outerposition',[0 0 1 1]);
histogram(log10(Dist_matrix(:,1)))
xlabel('Min distance to Aib (log10)')
ylabel('Count')
title(['A domain min distance distribution in the ',num2str(folders.num),' Trichoderma genomes (n=',num2str(sum(~isnan(Dist_matrix(:,1)))),')'])
xlim([-3,inf])
saveas(gcf,'./output/figure/A_min_dist_distribution_Trichoderma_log10.svg')
%% 
edge=0:0.005:0.15;
figure
histogram(Dist_matrix(Dist_matrix(:,1)<=0.15,1),edge)
hold on
xlabel('Min distance to Aib')
ylabel('Count')
loc_ylim=ylim;
line([0.075,0.075],loc_ylim,'Color','r','LineStyle','--')
hold off
saveas(gcf,'./output/figure/A_min_dist_distribution_Trichoderma_zoom_in.svg')
%% 
threshold=0.075;
% tabulate(folders.Genus(A_domain_folder_ids(Dist_matrix(:,1)<threshold)))
%% 
[uni_region_ids,~,ic]=unique(my_omains.region_ids(A_domain_index(Dist_matrix(:,1)<threshold)));
Aib_num=zeros(length(uni_region_ids),1);% the number of A domain which specify Aib
for i = 1:length(uni_region_ids)
    Aib_num(i)=sum(ic==i);
end
uni_sorted_product_type_str=my_regions.uni_sorted_product_type_str(uni_region_ids);
%% how many A domain in the region?
Aib_omain_ids=A_domain_index(Dist_matrix(:,1)<threshold);
locustag_Aib=my_omains.locustag(Aib_omain_ids);
NAD_domain=find(ismember(my_omains.domaintypelist,'NAD'));
TD_domain=find(ismember(my_omains.domaintypelist,'TD'));
PP_domain=find(ismember(my_omains.domaintypelist,'PP-binding'));
AT_domain=find(ismember(my_omains.domaintypelist,'AT'));
C_domain=find(ismember(my_omains.domaintypelist,'C'));
A_num=zeros(length(uni_region_ids),1);% the number of A domain in this region
A_singal_num=zeros(length(uni_region_ids),1);% the number of single (and only one domain) A domain in the single gene
C_num=zeros(length(uni_region_ids),1);% the number of c domain in this region
End_with_PP_NAD=zeros(length(uni_region_ids),1);% if contains TD domain, 1
End_with_PP_TD=zeros(length(uni_region_ids),1);% if ending with TD domain, 1
End_with_NAD=zeros(length(uni_region_ids),1);% if contains TD domain, 1
End_with_TD=zeros(length(uni_region_ids),1);% if ending with TD domain, 1
Aib_locustag=cell(length(uni_region_ids),1);% the locustag of A domain with Aib substrate 
Aib_locustag_A_ids=cell(length(uni_region_ids),1);% the locustag of A domain with Aib substrate | the index of A which specify Aib
uni_Aib_locustag=cell(length(uni_region_ids),1);% the unique locustag of A domain with Aib substrate 
BGC_architecture=cell(length(uni_region_ids),1);% architecture of BGC with Aib A domain
Start_with_Aib=zeros(length(uni_region_ids),1);% if there is one CDS which starts with A domain which specificity Aib, 1
AT_before_start_Aib=zeros(length(uni_region_ids),1); % if there is AT domain before the first A domain which specificity Aib, 1
Contain_one_AT=zeros(length(uni_region_ids),1); % if there is one and only one AT domain in this BGC
Module_iscomplete=zeros(length(uni_region_ids),1); % if #C~=#module, the module is incomplete, Module_iscomplete=1;
BGC_isfragmented=zeros(length(uni_region_ids),1); % if BGC containing modules is splited into many gene, BGC_isfragmented=1;
for i = 1:length(uni_region_ids)
    loc_typeid_mat=my_omains.typeid_mat(my_omains.region_ids==uni_region_ids(i),:);
    A_num(i)=sum(ismember(loc_typeid_mat,[A_domain,0],'rows'));
    loc_NAD_end=find(ismember(loc_typeid_mat,[NAD_domain,-1],'rows'));
    if ~isempty(loc_NAD_end)
        End_with_NAD(i)=1;
    end
    for j = 1:length(loc_NAD_end)
        if loc_typeid_mat(loc_NAD_end(j)-2,1)==PP_domain
            End_with_PP_NAD(i)=1;
        end
    end
    loc_TD_end=find(ismember(loc_typeid_mat,[TD_domain,-1],'rows'));
    if ~isempty(loc_TD_end)
        End_with_TD(i)=1;
    end
    for j = 1:length(loc_TD_end)
        if loc_typeid_mat(loc_TD_end(j)-2,1)==PP_domain
            End_with_PP_TD(i)=1;
        end
    end
    
    loc_locustag_Aib=my_omains.locustag(my_omains.region_ids==uni_region_ids(i));
    loc_omain_ids=find(my_omains.region_ids==uni_region_ids(i));
    loc_locustag_Aib=loc_locustag_Aib(loc_typeid_mat(:,2)==0);
    loc_omain_ids=loc_omain_ids(loc_typeid_mat(:,2)==0);
    loc_typeid_mat=loc_typeid_mat(loc_typeid_mat(:,2)==0,:);
    [uni_loc_locustag_Aib,~,ic1]=unique(loc_locustag_Aib,'stable');
    uni_loc_locustag_Aib_str=cell(length(uni_loc_locustag_Aib),1);
    if sum(loc_typeid_mat(:,1)==AT_domain)==1
        Contain_one_AT(i)=1;
    end
    loc_fragment=0;
    loc_uni_loc_locustag_Aib_index=1:length(uni_loc_locustag_Aib);
    remove_list=[];
    for j = 1:length(uni_loc_locustag_Aib)
        loc_index=find(ic1==j);
        loc_locustag_Aib_str=cell(length(loc_index),1);
        loc_typeid_mat_Aib=loc_typeid_mat(loc_index,1);
        if length(loc_typeid_mat_Aib)==1&&sum(loc_typeid_mat_Aib==A_domain)==1
            A_singal_num(i)=A_singal_num(i)+1;
            remove_list=[remove_list;j];
        elseif sum(loc_typeid_mat_Aib==A_domain)>0
            loc_fragment=loc_fragment+1;
        end
        C_num(i)=C_num(i)+sum(loc_typeid_mat_Aib==C_domain);
        for k = 1:length(loc_index)
            loc_locustag_Aib_str{k}=my_omains.domaintypelist{loc_typeid_mat(loc_index(k),1)};
            if ismember(loc_omain_ids(loc_index(k)),Aib_omain_ids)
                loc_locustag_Aib_str{k}=[loc_locustag_Aib_str{k},'(Aib)'];
            end
        end
        uni_loc_locustag_Aib_str(j)=join(loc_locustag_Aib_str,'+');
        uni_loc_locustag_Aib_str{j}=[uni_loc_locustag_Aib{j},'|',uni_loc_locustag_Aib_str{j}];
    end
    if loc_fragment>1
        BGC_isfragmented(i)=1;
    end
    BGC_architecture{i}=join(uni_loc_locustag_Aib_str,'; ');
    Aib_locustag(i)=join(locustag_Aib(ic==i),', ');
    loc_uni_locustag_Aib=unique(locustag_Aib(ic==i),'stable');
    uni_Aib_locustag(i)=join(loc_uni_locustag_Aib,', ');
    loc_Aib_index=Aib_omain_ids(ic==i);
    loc_Aib_locustag_A_ids=[];
    for j = 1:length(loc_uni_locustag_Aib)
        loc_A_index=find(my_omains.typeid_mat(:,1)==A_domain&my_omains.isdomain==1&ismember(my_omains.locustag,loc_uni_locustag_Aib{j})&my_omains.region_ids==uni_region_ids(i));
        [~,idx]=ismember(loc_Aib_index,loc_A_index);
        idx(idx==0)=[];
        for k = 1:length(idx)
            loc_Aib_locustag_A_ids=[loc_Aib_locustag_A_ids;{[loc_uni_locustag_Aib{j},'|',num2str(idx(k))]}];
            if idx(k)==1
                Start_with_Aib(i)=1;
            end
        end
        if Start_with_Aib(i)==1
            border=find(ismember(my_omains.locustag,loc_uni_locustag_Aib{j})&my_omains.region_ids==uni_region_ids(i),1);
            if any(my_omains.typeid_mat(border:loc_A_index(1)-1,1)==AT_domain)
                AT_before_start_Aib(i)=1;
            end
        end
    end
    assert(sum(ic==i)==length(loc_Aib_locustag_A_ids))
    Aib_locustag_A_ids(i)=join(loc_Aib_locustag_A_ids,', ');
    loc_uni_loc_locustag_Aib_index(remove_list)=[];
    if A_num(i)-A_singal_num(i)==C_num(i)
        if BGC_isfragmented(i)==1
            Module_iscomplete(i)=1;
        else
            loc_complete=zeros(length(loc_uni_loc_locustag_Aib_index),1);
            for j = 1:length(loc_uni_loc_locustag_Aib_index)
                loc_index=find(ic1==loc_uni_loc_locustag_Aib_index(j));
                loc_typeid_mat_Aib=loc_typeid_mat(loc_index,1);
                loc_C_index=find(loc_typeid_mat_Aib==C_domain);
                loc_A_index=find(loc_typeid_mat_Aib==A_domain);
                if loc_typeid_mat_Aib(1)~=A_domain&&loc_typeid_mat_Aib(end)~=C_domain&&all(loc_typeid_mat_Aib(loc_C_index+1)==A_domain)&&all(loc_typeid_mat_Aib(loc_A_index-1)==C_domain)
                    loc_complete(j)=1;
                end
            end
            if all(loc_complete==1)
                Module_iscomplete(i)=1;
            end
        end
    end
end
End_with_PP_TD_NAD=End_with_PP_NAD+End_with_PP_TD;
End_with_TD_NAD=End_with_NAD+End_with_TD;
Module_num=A_num-A_singal_num;
%% 
Known_BGC=cell(length(uni_region_ids),1);
real_Aib_num=cell(length(uni_region_ids),1);
Species=folders.Species(my_regions.folderid(uni_region_ids));
% Known_peptaibol=folders.Known_peptaibol(my_regions.folderid(uni_region_ids));
Known_peptaibol=cell(length(uni_region_ids),1);
foldername=my_regions.foldername(uni_region_ids);
filename=my_regions.filename(uni_region_ids);
for i = 1:length(foldername)
    loc_foldername=split(foldername{i},'/');
    foldername(i)=loc_foldername(end);
end
Species_short=cell(size(Species));
Strain=cell(size(Species));
for i = 1:length(Species)
    loc_str=split(Species{i},' ');
    Species_short{i}=['Trichoderma ',loc_str{2}];
    Strain{i}=strrep(Species{i},'T.','Trichoderma');
    Species{i}=loc_str{2};
end
%%
Aib_infor=table(Species_short,Aib_num,real_Aib_num,Module_num,Start_with_Aib,AT_before_start_Aib,Contain_one_AT,End_with_PP_TD,End_with_PP_NAD,End_with_PP_TD_NAD,End_with_TD,End_with_NAD,End_with_TD_NAD,uni_sorted_product_type_str,Known_peptaibol,Known_BGC,BGC_architecture,Aib_locustag_A_ids,uni_Aib_locustag,foldername,filename,Species,Strain,Module_iscomplete,BGC_isfragmented);
writetable(Aib_infor,'Aib_infor_T.xlsx','Sheet','infor')
%% 
Peptaibol_check=readmatrix('Aib_infor_T.xlsx','Sheet','infor0','Range','Z2:Z321');
%% after checking known peptaibol in infor sheet
peptaibol_check=readcell('Aib_infor_T.xlsx','Sheet','infor0','Range','O:O');
genome_check=readcell('Aib_infor_T.xlsx','Sheet','infor0','Range','T:T');
peptaibol_check(1)=[];
genome_check(1)=[];
%% 
Known_peptaibol=peptaibol_check;
Known_peptaibol_index=ones(length(Known_peptaibol),1);
for i = 1:length(Known_peptaibol)
    if ismissing(Known_peptaibol{i})
        Known_peptaibol{i}=[];
        Known_peptaibol_index(i)=0;
    end
end
%% 
Known_peptaibol_uni_region_ids=uni_region_ids(Known_peptaibol_index==1);
Known_peptaibol_struct=[];
Known_peptaibol_struct.num=length(Known_peptaibol_uni_region_ids);
Known_peptaibol_struct.peptaibol_index=find(Known_peptaibol_index==1);
Known_peptaibol_struct.peptaibol_check=Peptaibol_check(Known_peptaibol_index==1);
Known_peptaibol_struct.Module_iscomplete=Module_iscomplete(Known_peptaibol_index==1);
Known_peptaibol_struct.BGC_isfragmented=BGC_isfragmented(Known_peptaibol_index==1);
Known_peptaibol_struct.Module_num=Module_num(Known_peptaibol_index==1);
Known_peptaibol_struct.peptaibol_name=Known_peptaibol(Known_peptaibol_index==1);
Known_peptaibol_struct.uni_Aib_locustag=uni_Aib_locustag(Known_peptaibol_index==1);
Known_peptaibol_struct.foldername=foldername(Known_peptaibol_index==1);
Known_peptaibol_struct.filename=filename(Known_peptaibol_index==1);
Known_peptaibol_struct.Species=Species(Known_peptaibol_index==1);
Known_peptaibol_struct.Strain=Strain(Known_peptaibol_index==1);
%% 
for i = 1:Known_peptaibol_struct.num
    Known_peptaibol_struct.original_ntseq(i,1)=join(my_omains.seq_ntaa(ismember(my_omains.locustag,Known_peptaibol_struct.uni_Aib_locustag{i})&(my_omains.region_ids==Known_peptaibol_uni_region_ids(i)),1),'');
    Known_peptaibol_struct.NT_seqs{i,1}=Known_peptaibol_struct.original_ntseq{i};
    Known_peptaibol_struct.AA_seqs{i,1}=nt2aa(Known_peptaibol_struct.original_ntseq{i},'AlternativeStartCodons',false,'ACGTOnly',false);
    X_AA_index=strfind(Known_peptaibol_struct.AA_seqs{i},'X');
    X_NT_index=[3*X_AA_index-2,3*X_AA_index-1,3*X_AA_index];
    Known_peptaibol_struct.AA_seqs{i}(X_AA_index)=[];
    Known_peptaibol_struct.NT_seqs{i}(X_NT_index)=[];
    Known_peptaibol_struct.AA_seq_len(i,1)=length(Known_peptaibol_struct.AA_seqs{i});
end
%% 
Known_peptaibol_struct.num=Known_peptaibol_struct.num+1;
Known_peptaibol_struct.Module_num=[Known_peptaibol_struct.Module_num;18];
Known_peptaibol_struct.peptaibol_name=[Known_peptaibol_struct.peptaibol_name;{'PlePS/NPS1tp'}];
Known_peptaibol_struct.uni_Aib_locustag=[Known_peptaibol_struct.uni_Aib_locustag;{'AQV12033.1'}];
Known_peptaibol_struct.foldername=[Known_peptaibol_struct.foldername;{''}];
Known_peptaibol_struct.filename=[Known_peptaibol_struct.filename;{'PlePS_T. pleuroti'}];
Known_peptaibol_struct.Species=[Known_peptaibol_struct.Species;{'pleuroti'}];
Known_peptaibol_struct.Strain=[Known_peptaibol_struct.Strain;{'Trichoderma pleuroti TPhu1'}];
Known_peptaibol_struct.peptaibol_check=[Known_peptaibol_struct.peptaibol_check;1];
Known_peptaibol_struct.Module_iscomplete=[Known_peptaibol_struct.Module_iscomplete;1];
Known_peptaibol_struct.BGC_isfragmented=[Known_peptaibol_struct.BGC_isfragmented;0];
%% 
NT_files={'PlePS_T. pleuroti.fasta'};
NT_seqs=cell(length(NT_files),1);
NT_AA_seqs=cell(length(NT_files),1);
for i = 1:length(NT_files)
    NT_AA_seqs{i}=fastaread(NT_files{i});
    NT_seqs{i}=fastaread(['cds-',strrep(NT_files{i},'.fasta','.txt')]);
    Known_peptaibol_struct.original_ntseq=[Known_peptaibol_struct.original_ntseq;{NT_seqs{i}.Sequence}];
    % 
end
% remove stop codon
for i = 1:length(NT_seqs)
    if length(NT_seqs{i}.Sequence)/3-1==length(NT_AA_seqs{i}.Sequence)
        NT_seqs{i}.Sequence(end-2:end)=[];
    end
end
for i = 1:length(NT_seqs)
    assert(length(NT_seqs{i}.Sequence)/3==length(NT_AA_seqs{i}.Sequence))
end
% remove X in the aa sequence because X is removed in the NRPS motif Finder
% result
for i = 1:length(NT_seqs)
    X_index=strfind(NT_AA_seqs{i}.Sequence,'X');
    X_AA_index=[3*X_index-2,3*X_index-1,3*X_index];
    NT_AA_seqs{i}.Sequence(X_index)=[];
    NT_seqs{i}.Sequence(X_AA_index)=[];
    assert(length(NT_seqs{i}.Sequence)/3==length(NT_AA_seqs{i}.Sequence))
    Known_peptaibol_struct.NT_seqs=[Known_peptaibol_struct.NT_seqs;{NT_seqs{i}.Sequence}];
    Known_peptaibol_struct.AA_seqs=[Known_peptaibol_struct.AA_seqs;{NT_AA_seqs{i}.Sequence}];
    Known_peptaibol_struct.AA_seq_len=[Known_peptaibol_struct.AA_seq_len;length(NT_AA_seqs{i}.Sequence)];
end
%% 
for i = 1:Known_peptaibol_struct.num
    loc_seq=[];
    loc_seq.Header=num2str(i);
    loc_seq.Sequence=Known_peptaibol_struct.AA_seqs{i};
    if i == 56
        Known_peptaibol_struct.result{i,1} = Find_NRPS_motif_module_pfam_HRL(loc_seq,[],[],{'Aalpha','G','Talpha'},0,0.5,[]);% For one Tex2, use smaller threshold
    else
        Known_peptaibol_struct.result{i,1} = Find_NRPS_motif_module_pfam_HRL(loc_seq,[],[],{'Aalpha','G','Talpha'},0,[],[]);
    end
end
%% 
for i = 1:Known_peptaibol_struct.num
    Known_peptaibol_struct.Pred_module_num(i,1)=sum(Known_peptaibol_struct.result{i}.domain_list==2);
end
assert(all(Known_peptaibol_struct.Pred_module_num==Known_peptaibol_struct.Module_num));
%% 
Known_peptaibol_struct.seq_matrix=[];
Known_peptaibol_struct.species_list=[];
Known_peptaibol_struct.module_list=[];
for i = 1:length(Known_peptaibol_struct.result)
    loc_index=find(ismember(Known_peptaibol_struct.result{i}.motifid_mat,[2,1],'rows'));
    for j = 1:length(loc_index)
        Known_peptaibol_struct.species_list=[Known_peptaibol_struct.species_list;i];
        Known_peptaibol_struct.module_list=[Known_peptaibol_struct.module_list;j];
        Known_peptaibol_struct.seq_matrix=[Known_peptaibol_struct.seq_matrix;[Known_peptaibol_struct.result{i}.seq_list(loc_index(j):loc_index(j)+23)]'];
    end
end
Known_peptaibol_struct.motifid_mat=Known_peptaibol_struct.result{i}.motifid_mat(loc_index(j):loc_index(j)+23,:);
%% 
NT_seq_matrix=cell(size(Known_peptaibol_struct.seq_matrix));
for i = 1:length(Known_peptaibol_struct.NT_seqs)
    loc_index=find(Known_peptaibol_struct.species_list==i);
    for j = 1:length(loc_index)
        loc_A_seq=join(Known_peptaibol_struct.seq_matrix(loc_index(j),:),'');
        loc_A_index=strfind(Known_peptaibol_struct.AA_seqs{i},loc_A_seq{1});
        if ~isempty(loc_A_index)
            motif_len=cellfun(@length,Known_peptaibol_struct.seq_matrix(loc_index(j),:));
            for k = 1:length(motif_len)
                % NT_seq_matrix{loc_index(j),k}=NT_AA_seqs{i}.Sequence(loc_A_index+sum(motif_len(1:k-1)):loc_A_index+sum(motif_len(1:k))-1);
                NT_seq_matrix{loc_index(j),k}=upper(Known_peptaibol_struct.NT_seqs{i}((loc_A_index+sum(motif_len(1:k-1)))*3-2:(loc_A_index+sum(motif_len(1:k))-1)*3));
                loc_NT2AA=nt2aa(NT_seq_matrix{loc_index(j),k},'AlternativeStartCodons',false,'ACGTOnly',false);
                if isempty(Known_peptaibol_struct.seq_matrix{loc_index(j),k})&&isempty(loc_NT2AA)
                    loc_NT2AA=Known_peptaibol_struct.seq_matrix{loc_index(j),k};
                end
                if ~strcmp(Known_peptaibol_struct.seq_matrix{loc_index(j),k},loc_NT2AA)
                    fprintf('The AA is %s, the NT is %s, the NT2AA is %s\n',Known_peptaibol_struct.seq_matrix{loc_index(j),k},NT_seq_matrix{loc_index(j),k},loc_NT2AA)
                end
            end
        else
            j % the 17th module of 73 peptaibol has problem,but 73 is dead. Known_peptaibol_struct.dead_peptaibol(73)==1. omit this problem
        end
    end
end
Known_peptaibol_struct.NT_seq_matrix=NT_seq_matrix;
%% load substrate specifity
substrate_raw=readcell('Peptaibol NRPS and Peptaibols.xlsx','Sheet','PSs and corresponding peptaibol');

for i = 1:size(substrate_raw,1)
    for j = 1:size(substrate_raw,2)
        if ismissing(substrate_raw{i,j})
            substrate_raw{i,j}=[];
        end
    end
end

peptaibol_name_index=[];
for i = 1:size(substrate_raw,1)
    if strcmp(substrate_raw{i,3},'Peptaibol synthetase')
        peptaibol_name_index=[peptaibol_name_index;i];
    end
end
peptaibol_name_list=substrate_raw(peptaibol_name_index+1,3);
peptaibol_num_list=cell2mat(substrate_raw(peptaibol_name_index+1,2));
peptaibol_best_strain_list=substrate_raw(peptaibol_name_index+1,1);
%% 
known_species_substrate_specificity=cell(length(peptaibol_num_list),1);
for i = 1:length(peptaibol_num_list)
    if i~=length(peptaibol_num_list)
        if strcmp(peptaibol_name_list{i},'AtvPS/Pbs1')
            known_species_substrate_specificity{i}=substrate_raw(peptaibol_name_index(i)+3:peptaibol_name_index(i+1)-1,5:5+peptaibol_num_list(i));
            known_species_substrate_specificity{i}(:,6)=[];% skip one unknown module
        elseif strcmp(peptaibol_name_list{i},'ArdPS1')%ArdPS1 only has one product, but more than 3 species produce it.
            known_species_substrate_specificity{i}=substrate_raw(peptaibol_name_index(i)+3,5:4+peptaibol_num_list(i));
        else
            known_species_substrate_specificity{i}=substrate_raw(peptaibol_name_index(i)+3:peptaibol_name_index(i+1)-1,5:4+peptaibol_num_list(i));
        end
    else
        known_species_substrate_specificity{i}=substrate_raw(peptaibol_name_index(i)+3:end,5:4+peptaibol_num_list(i));
    end
end
%% 
peptaibol_name_list=[peptaibol_name_list;{'HypoPS1'}];
peptaibol_num_list=[peptaibol_num_list;19];
known_species_substrate_specificity=[known_species_substrate_specificity;{readcell('trichohypolin structure.xlsx','Sheet','trichohypolins 3 letters','Range','B3:T46')}];
peptaibol_best_strain_list=[peptaibol_best_strain_list;{[]};{'Trichoderma hypoxylon'}];
Known_peptaibol_struct.peptaibol_name_list=peptaibol_name_list;
Known_peptaibol_struct.peptaibol_num_list=peptaibol_num_list;
Known_peptaibol_struct.peptaibol_best_strain_list=peptaibol_best_strain_list;
%% 
for i = 1:length(known_species_substrate_specificity)
    for j = 1:size(known_species_substrate_specificity{i},1)
        for k = 1:size(known_species_substrate_specificity{i},2)
            if strcmp(known_species_substrate_specificity{i}{j,k},'Ile')||strcmp(known_species_substrate_specificity{i}{j,k},'Ileol')||strcmp(known_species_substrate_specificity{i}{j,k},'Leu')||strcmp(known_species_substrate_specificity{i}{j,k},'Leuol')||strcmp(known_species_substrate_specificity{i}{j,k},'Lxxol')
                known_species_substrate_specificity{i}{j,k}='Lxx';
            elseif strcmp(known_species_substrate_specificity{i}{j,k},'Val')||strcmp(known_species_substrate_specificity{i}{j,k},'Valol')||strcmp(known_species_substrate_specificity{i}{j,k},'Vxxol')||strcmp(known_species_substrate_specificity{i}{j,k},'Iva')
                known_species_substrate_specificity{i}{j,k}='Vxx';
            elseif strcmp(known_species_substrate_specificity{i}{j,k},'Glu(Me)')
                known_species_substrate_specificity{i}{j,k}='Glu(OMe)';
            elseif strcmp(known_species_substrate_specificity{i}{j,k},'Pheol')
                known_species_substrate_specificity{i}{j,k}='Phe';
            elseif strcmp(known_species_substrate_specificity{i}{j,k},'Trpol')
                known_species_substrate_specificity{i}{j,k}='Trp';
            end
        end
    end
end
%% 
all_substrate_list=[];
for j = 1:length(known_species_substrate_specificity)
    for i =1:size(known_species_substrate_specificity{j},2)
        all_substrate_list=[all_substrate_list;known_species_substrate_specificity{j}(:,i)];
    end
end
unique_substrate=unique(all_substrate_list);
for i = 1:length(unique_substrate)
    if strcmp(unique_substrate(i),'Empty')
        unique_substrate(i)=[];
        break
    end
end
unique_substrate=[unique_substrate;{'Empty'}];
%% 
known_species_substrate_specificity_n=cell(length(known_species_substrate_specificity),1);
known_species_substrate_specificity_freq=cell(length(known_species_substrate_specificity),1);
for k = 1:length(known_species_substrate_specificity)
    substrate_n_matrix=[];
    for i = 1:size(known_species_substrate_specificity{k},2)
        loc_col=known_species_substrate_specificity{k}(:,i);
        loc_num=zeros(size(known_species_substrate_specificity{k},1),1);
        for j = 1:length(unique_substrate)
            loc_num=loc_num+j*strcmp(loc_col,unique_substrate(j));
        end
        substrate_n_matrix=[substrate_n_matrix,loc_num];
    end
    known_species_substrate_specificity_n{k}=substrate_n_matrix;
    substrate_tab=zeros(length(unique_substrate),size(known_species_substrate_specificity{k},2));
    col_n=size(known_species_substrate_specificity{k},1);
    for i = 1:size(known_species_substrate_specificity{k},2)
        for j = 1:length(unique_substrate)
            substrate_tab(j,i)=sum(substrate_n_matrix(:,i)==j)/col_n;
        end
    end
    known_species_substrate_specificity_freq{k}=substrate_tab';
end
%% 
specificity_freq=[];
for i = 1:Known_peptaibol_struct.num
    specificity_freq=[specificity_freq;known_species_substrate_specificity_freq{ismember(peptaibol_name_list,Known_peptaibol_struct.peptaibol_name{i})}];
end
substrate_mode=zeros(size(specificity_freq,1),1);
for i = 1:size(specificity_freq,1)
    [~,substrate_mode(i)]=max(specificity_freq(i,:));
end
Known_peptaibol_struct.specificity_freq=specificity_freq;
Known_peptaibol_struct.substrate_mode=substrate_mode;
%% 
Known_peptaibol_struct.seq_len_matrix=cellfun(@length,Known_peptaibol_struct.seq_matrix);
Known_peptaibol_struct.unique_substrate=unique_substrate;
%% unknown BGC but known substrates
substrate_raw=readcell('Peptaibol NRPS and Peptaibols.xlsx','Sheet','Peptaibols but no NRPS proposed');
%%
for i = 1:size(substrate_raw,1)
    for j = 1:size(substrate_raw,2)
        if ismissing(substrate_raw{i,j})
            substrate_raw{i,j}=[];
        end
    end
end
%% 
peptaibol_name_index=[];
for i = 1:size(substrate_raw,1)
    if strcmp(substrate_raw{i,3},'Compound Name')
        peptaibol_name_index=[peptaibol_name_index;i];
    end
end
peptaibol_num_list=cell2mat(substrate_raw(peptaibol_name_index+1,1));
%% 
unknown_species_substrate_specificity=cell(length(peptaibol_num_list),1);
for i = 1:length(peptaibol_num_list)
    if i~=length(peptaibol_num_list)
        unknown_species_substrate_specificity{i}=substrate_raw(peptaibol_name_index(i)+1:peptaibol_name_index(i+1)-1,5:4+peptaibol_num_list(i));
    else
        unknown_species_substrate_specificity{i}=substrate_raw(peptaibol_name_index(i)+1:end,5:4+peptaibol_num_list(i));
    end
end
%% 
for i = 1:length(unknown_species_substrate_specificity)
    for j = 1:size(unknown_species_substrate_specificity{i},1)
        for k = 1:size(unknown_species_substrate_specificity{i},2)
            if strcmp(unknown_species_substrate_specificity{i}{j,k},'Ile')||strcmp(unknown_species_substrate_specificity{i}{j,k},'Ileol')||strcmp(unknown_species_substrate_specificity{i}{j,k},'Leu')||strcmp(unknown_species_substrate_specificity{i}{j,k},'Leuol')||strcmp(unknown_species_substrate_specificity{i}{j,k},'Lxxol')
                unknown_species_substrate_specificity{i}{j,k}='Lxx';
            elseif strcmp(unknown_species_substrate_specificity{i}{j,k},'Val')||strcmp(unknown_species_substrate_specificity{i}{j,k},'Valol')||strcmp(unknown_species_substrate_specificity{i}{j,k},'Vxxol')||strcmp(unknown_species_substrate_specificity{i}{j,k},'Iva')
                unknown_species_substrate_specificity{i}{j,k}='Vxx';
            elseif strcmp(unknown_species_substrate_specificity{i}{j,k},'Glu(Me)')
                unknown_species_substrate_specificity{i}{j,k}='Glu(OMe)';
            elseif strcmp(unknown_species_substrate_specificity{i}{j,k},'Pheol')
                unknown_species_substrate_specificity{i}{j,k}='Phe';
            elseif strcmp(unknown_species_substrate_specificity{i}{j,k},'Trpol')
                unknown_species_substrate_specificity{i}{j,k}='Trp';
            end
        end
    end
end
%% 
for j = 1:length(unknown_species_substrate_specificity)
    for i =1:size(unknown_species_substrate_specificity{j},2)
        all_substrate_list=[all_substrate_list;unknown_species_substrate_specificity{j}(:,i)];
    end
end
unique_substrate_unknown=unique(all_substrate_list);
%% 
unique_substrate=[Known_peptaibol_struct.unique_substrate;{'Lys'}];
assert(all(ismember(unique_substrate_unknown,unique_substrate)))
%% 
unknown_species_substrate_specificity_n=cell(length(unknown_species_substrate_specificity),1);
unknown_species_substrate_specificity_freq=cell(length(unknown_species_substrate_specificity),1);
for k = 1:length(unknown_species_substrate_specificity)
    substrate_n_matrix=[];
    for i = 1:size(unknown_species_substrate_specificity{k},2)
        loc_col=unknown_species_substrate_specificity{k}(:,i);
        loc_num=zeros(size(unknown_species_substrate_specificity{k},1),1);
        for j = 1:length(unique_substrate)
            loc_num=loc_num+j*strcmp(loc_col,unique_substrate(j));
        end
        substrate_n_matrix=[substrate_n_matrix,loc_num];
    end
    unknown_species_substrate_specificity_n{k}=substrate_n_matrix;
    substrate_tab=zeros(length(unique_substrate),size(unknown_species_substrate_specificity{k},2));
    col_n=size(unknown_species_substrate_specificity{k},1);
    for i = 1:size(unknown_species_substrate_specificity{k},2)
        for j = 1:length(unique_substrate)
            substrate_tab(j,i)=sum(substrate_n_matrix(:,i)==j)/col_n;
        end
    end
    unknown_species_substrate_specificity_freq{k}=substrate_tab';
end
%% 
unique_substrate=Known_peptaibol_struct.unique_substrate;
specificity_freq=[];
for i = 1:length(known_species_substrate_specificity_freq)
    specificity_freq=[specificity_freq;known_species_substrate_specificity_freq{i}];
end
loc_path='output/figure/substrate_promiscuity/18peptaibol/';
if 1
    unique_substrate=[unique_substrate;{'Lys'}];
    specificity_freq=[specificity_freq,zeros(size(specificity_freq,1),1)];
    for i = 1:length(unknown_species_substrate_specificity_freq)
        specificity_freq=[specificity_freq;unknown_species_substrate_specificity_freq{i}];
    end
    loc_path='output/figure/substrate_promiscuity/26peptaibol/';
end
node_size=sum(specificity_freq)';
node_infor=readcell([loc_path,'node.xlsx'],'Sheet','Sheet1','Range','A:C');
node=[node_infor,[{'Size'};num2cell(node_size)]];
 
writecell(node,[loc_path,'node.xlsx'],'Sheet','Sheet1');

edge=[];
for i = 1:length(unique_substrate)
    for j = i+1:length(unique_substrate)
        weight=0;
        for k = 1:size(specificity_freq,1)
            if specificity_freq(k,i)*specificity_freq(k,j)>0
                weight=weight+(specificity_freq(k,i)+specificity_freq(k,j))/2;
            end
        end
        if weight>0
            edge=[edge;[i,j,weight]];
        end
    end
end
edge_name=[unique_substrate(edge(:,1)),unique_substrate(edge(:,2))];
edge_name=[edge_name,repmat({'FLASE'},size(edge,1),1),num2cell(edge(:,3))];
edge_name=[{'source','target','directed','weight'};edge_name];
writecell(edge_name,[loc_path,'edge.xlsx']);
%% 
known_species_substrate_A_n=0;
for i = 1:length(known_species_substrate_specificity)
    known_species_substrate_A_n=known_species_substrate_A_n+size(known_species_substrate_specificity{i},2);
end
unknown_species_substrate_A_n=0;
for i = 1:length(unknown_species_substrate_specificity)
    unknown_species_substrate_A_n=unknown_species_substrate_A_n+size(unknown_species_substrate_specificity{i},2);
end
%% 
seq_len_mode=mode(Known_peptaibol_struct.seq_len_matrix);
%% 
motif_index=[3:2:11,15:2:23];%A1-A10
dead_A=zeros(length(Known_peptaibol_struct.species_list),1);
for i = 1:length(Known_peptaibol_struct.species_list)
    if ~all(Known_peptaibol_struct.seq_len_matrix(i,motif_index)==seq_len_mode(motif_index))
        dead_A(i)=1;
    end
end
dead_peptaibol=unique(Known_peptaibol_struct.species_list(dead_A==1));
Known_peptaibol_struct.dead_peptaibol=zeros(length(Known_peptaibol_struct.Species),1);
Known_peptaibol_struct.dead_peptaibol(dead_peptaibol)=1;
Known_peptaibol_struct.dead_peptaibol(end)=0;%PlePS/NPS1tp only has one species, and shorter in A6, A7, longer in A8, A3-A6 is ok(one gap in A6)
%% check A domain dist within one type PS
good_peptaibol=[];
good_peptaibol.num=length(Known_peptaibol_struct.peptaibol_name_list);
for i = 1:good_peptaibol.num
    good_peptaibol.peptaibol_name{i,1}=Known_peptaibol_struct.peptaibol_name_list{i};
    good_peptaibol.species_index{i,1}=find(ismember(Known_peptaibol_struct.peptaibol_name,Known_peptaibol_struct.peptaibol_name_list{i})&Known_peptaibol_struct.dead_peptaibol==0&Known_peptaibol_struct.Module_iscomplete==1&Known_peptaibol_struct.BGC_isfragmented==0);
    if ~isempty(good_peptaibol.species_index{i,1})
        good_peptaibol.A_num(i,1)=Known_peptaibol_struct.peptaibol_num_list(i);
        for j = 1:good_peptaibol.A_num(i)
            for k = 7:11 %A3-A6
                good_peptaibol.seq_matrix{i,1}{j,k-6}=Known_peptaibol_struct.seq_matrix(ismember(Known_peptaibol_struct.species_list,good_peptaibol.species_index{i})&Known_peptaibol_struct.module_list==j,k);
                good_peptaibol.NT_seq_matrix{i,1}{j,k-6}=Known_peptaibol_struct.NT_seq_matrix(ismember(Known_peptaibol_struct.species_list,good_peptaibol.species_index{i})&Known_peptaibol_struct.module_list==j,k);
            end
            %merge G motif due to some sequences don't have G-motif
            loc_seq=vertcat(Known_peptaibol_struct.seq_matrix(ismember(Known_peptaibol_struct.species_list,good_peptaibol.species_index{i})&Known_peptaibol_struct.module_list==j,12:14));
            merge_seq=cell(size(loc_seq,1),1);
            for k = 1:size(loc_seq,1)
                merge_seq{k}=[loc_seq{k,1} loc_seq{k,2} loc_seq{k,3}];
            end
            good_peptaibol.seq_matrix{i,1}{j,6}=merge_seq;
            loc_seq=vertcat(Known_peptaibol_struct.NT_seq_matrix(ismember(Known_peptaibol_struct.species_list,good_peptaibol.species_index{i})&Known_peptaibol_struct.module_list==j,12:14));
            merge_seq=cell(size(loc_seq,1),1);
            for k = 1:size(loc_seq,1)
                merge_seq{k}=[loc_seq{k,1} loc_seq{k,2} loc_seq{k,3}];
            end
            good_peptaibol.NT_seq_matrix{i,1}{j,6}=merge_seq;
            for k = 15 %A3-A6
                good_peptaibol.seq_matrix{i,1}{j,k-8}=Known_peptaibol_struct.seq_matrix(ismember(Known_peptaibol_struct.species_list,good_peptaibol.species_index{i})&Known_peptaibol_struct.module_list==j,k);
                good_peptaibol.NT_seq_matrix{i,1}{j,k-8}=Known_peptaibol_struct.NT_seq_matrix(ismember(Known_peptaibol_struct.species_list,good_peptaibol.species_index{i})&Known_peptaibol_struct.module_list==j,k);
            end
        end
    else
        i
    end
end
%% 
good_peptaibol.msa_matrix=[];
loc_pwd=pwd;
if exist('tmp_24601','dir')==7
    rmdir('tmp_24601', 's')
end
mkdir tmp_24601 % bulid a temporary folder
for i = 1:good_peptaibol.num
    for j = 1:size(good_peptaibol.seq_matrix{i},1)
        for k = 1:size(good_peptaibol.seq_matrix{i},2)
            if mod(k,2)~=0%motif
                good_peptaibol.msa_matrix{i,1}{j,k}=[];
                for kk = 1:length(good_peptaibol.seq_matrix{i,1}{j,k})
                    good_peptaibol.msa_matrix{i,1}{j,k}=[good_peptaibol.msa_matrix{i,1}{j,k};good_peptaibol.seq_matrix{i,1}{j,k}{kk}];
                end
            else%intermotif
                if length(good_peptaibol.seq_matrix{i,1}{j,k})>1
                    test_header = cellstr(num2str([1:length(good_peptaibol.seq_matrix{i,1}{j,k})]'));
                    fastawrite([pwd,'/tmp_24601/',num2str(i),'_',num2str(j),'_',num2str(k),'.fasta'],test_header,good_peptaibol.seq_matrix{i,1}{j,k});
                    command = ['clustalo -i ',loc_pwd,'/tmp_24601/',num2str(i),'_',num2str(j),'_',num2str(k),'.fasta -o ',loc_pwd,'/tmp_24601/',num2str(i),'_',num2str(j),'_',num2str(k),'_1.fasta --outfmt fasta --output-order input-order'];
                    system(command);
                    good_peptaibol.msa_matrix{i,1}{j,k}=[];
                    raw = fastaread([pwd,'/tmp_24601/',num2str(i),'_',num2str(j),'_',num2str(k),'_1.fasta']);
                    for ii = 1:length(raw)
                        good_peptaibol.msa_matrix{i,1}{j,k}=[good_peptaibol.msa_matrix{i,1}{j,k};raw(ii).Sequence];
                    end
                else
                    good_peptaibol.msa_matrix{i,1}{j,k}=good_peptaibol.seq_matrix{i,1}{j,k}{1};
                end
            end
        end
    end
end
rmdir('tmp_24601', 's')
%% 
good_peptaibol.dist_matrix=cell(size(good_peptaibol.msa_matrix));% A3-A6 intermotif motif split
for i = 1:good_peptaibol.num
    good_peptaibol.dist_matrix{i}=cell(size(good_peptaibol.msa_matrix{i}));
    for j = 1:size(good_peptaibol.msa_matrix{i},1)
        for k = 1:size(good_peptaibol.msa_matrix{i},2)
            if size(good_peptaibol.msa_matrix{i}{j,k},1)>1
                good_peptaibol.dist_matrix{i}{j,k}=seqpdist(good_peptaibol.msa_matrix{i}{j,k},'ScoringMatrix','BLOSUM62','Method','alignment-score','SquareForm',true);
            else
                good_peptaibol.dist_matrix{i}{j,k}=0;
            end
        end
    end
end
%% 
good_peptaibol.dist_list=cell(size(good_peptaibol.msa_matrix)); % A3-A6 the whole
good_peptaibol.p_dist_list=cell(size(good_peptaibol.msa_matrix)); 
good_peptaibol.dist_max_list=cell(size(good_peptaibol.msa_matrix));
good_peptaibol.p_dist_max_list=cell(size(good_peptaibol.msa_matrix));
for i = 1:good_peptaibol.num
    good_peptaibol.dist_list{i}=cell(size(good_peptaibol.msa_matrix{i},1),1);
    good_peptaibol.p_dist_list{i}=cell(size(good_peptaibol.msa_matrix{i},1),1);
    good_peptaibol.dist_max_list{i}=zeros(size(good_peptaibol.msa_matrix{i},1),1);
    good_peptaibol.p_dist_max_list{i}=zeros(size(good_peptaibol.msa_matrix{i},1),1);
    for j = 1:size(good_peptaibol.msa_matrix{i},1)
        if size(good_peptaibol.msa_matrix{i}{1,1},1)==1
            good_peptaibol.dist_list{i}{j,1}=0;
            good_peptaibol.p_dist_list{i}{j,1}=0;
        else
            loc_seq=[];
            for k = 1:size(good_peptaibol.msa_matrix{i},2)
                loc_seq=[loc_seq,good_peptaibol.msa_matrix{i}{j,k}];
            end
            good_peptaibol.dist_list{i}{j,1}=seqpdist(loc_seq,'ScoringMatrix','BLOSUM62','Method','alignment-score','SquareForm',true);
            good_peptaibol.p_dist_list{i}{j,1}=seqpdist(loc_seq,'ScoringMatrix','BLOSUM62','Method','p-distance','SquareForm',true);
        end
        good_peptaibol.dist_max_list{i}(j)=max(good_peptaibol.dist_list{i}{j},[],'all');
        good_peptaibol.p_dist_max_list{i}(j)=max(good_peptaibol.p_dist_list{i}{j},[],'all');
    end
end
%% 
good_peptaibol.species_best_index=zeros(good_peptaibol.num,1);
for i = 1:good_peptaibol.num
    if ~isempty(Known_peptaibol_struct.peptaibol_best_strain_list{i})
        good_peptaibol.species_best_index(i,1)=good_peptaibol.species_index{i}(ismember(Known_peptaibol_struct.Strain(good_peptaibol.species_index{i}),Known_peptaibol_struct.peptaibol_best_strain_list{i}));
    else
        if length(good_peptaibol.species_index{i})==1
            good_peptaibol.species_best_index(i,1)=good_peptaibol.species_index{i};
        else
            loc_dist_sum=0;
            for j = 1:length(good_peptaibol.p_dist_list{i})
                loc_dist_sum=loc_dist_sum+sum(good_peptaibol.p_dist_list{i}{j});
            end
            [~,I]=min(loc_dist_sum);
            good_peptaibol.species_best_index(i,1)=good_peptaibol.species_index{i}(I);
        end
    end
end
%% 
best_peptaibol=[];
best_peptaibol.num=length(Known_peptaibol_struct.peptaibol_name_list);
best_peptaibol.species_best_index=good_peptaibol.species_best_index;
best_peptaibol.Strain_used=Known_peptaibol_struct.Strain(good_peptaibol.species_best_index);
best_peptaibol.A_seqs_list=cell(best_peptaibol.num,1);
for i = 1:best_peptaibol.num
    best_peptaibol.peptaibol_name{i,1}=Known_peptaibol_struct.peptaibol_name_list{i};
    best_peptaibol.A_num(i,1)=Known_peptaibol_struct.peptaibol_num_list(i);
    best_peptaibol.A_seqs_list{i}=join(Known_peptaibol_struct.seq_matrix(ismember(Known_peptaibol_struct.species_list,best_peptaibol.species_best_index(i)),:),'');
    for j = 1:best_peptaibol.A_num(i)
        for k = 7:11 %A3-A6
            best_peptaibol.seq_matrix{i,1}{j,k-6}=Known_peptaibol_struct.seq_matrix{ismember(Known_peptaibol_struct.species_list,best_peptaibol.species_best_index(i))&Known_peptaibol_struct.module_list==j,k};
            best_peptaibol.NT_seq_matrix{i,1}{j,k-6}=Known_peptaibol_struct.NT_seq_matrix{ismember(Known_peptaibol_struct.species_list,best_peptaibol.species_best_index(i))&Known_peptaibol_struct.module_list==j,k};
        end
        %merge G motif due to some sequences don't have G-motif
        loc_seq=Known_peptaibol_struct.seq_matrix(ismember(Known_peptaibol_struct.species_list,best_peptaibol.species_best_index(i))&Known_peptaibol_struct.module_list==j,12:14);
        merge_seq=[];
        for k = 1:size(loc_seq,2)
            merge_seq=[merge_seq,loc_seq{k}];
        end
        best_peptaibol.seq_matrix{i,1}{j,6}=merge_seq;
        loc_seq=Known_peptaibol_struct.NT_seq_matrix(ismember(Known_peptaibol_struct.species_list,best_peptaibol.species_best_index(i))&Known_peptaibol_struct.module_list==j,12:14);
        merge_seq=[];
        for k = 1:size(loc_seq,2)
            merge_seq=[merge_seq,loc_seq{k}];
        end
        best_peptaibol.NT_seq_matrix{i,1}{j,6}=merge_seq;
        for k = 15 %A3-A6
            best_peptaibol.seq_matrix{i,1}{j,k-8}=Known_peptaibol_struct.seq_matrix{ismember(Known_peptaibol_struct.species_list,best_peptaibol.species_best_index(i))&Known_peptaibol_struct.module_list==j,k};
            best_peptaibol.NT_seq_matrix{i,1}{j,k-8}=Known_peptaibol_struct.NT_seq_matrix{ismember(Known_peptaibol_struct.species_list,best_peptaibol.species_best_index(i))&Known_peptaibol_struct.module_list==j,k};
        end
    end
end
best_peptaibol.specificity_freq=known_species_substrate_specificity_freq;
best_peptaibol.specificity_freq_matrix=[];
for i = 1:best_peptaibol.num
    best_peptaibol.specificity_freq_matrix=[best_peptaibol.specificity_freq_matrix;best_peptaibol.specificity_freq{i}];
end
substrate_mode=zeros(size(best_peptaibol.specificity_freq_matrix,1),1);
for i = 1:size(best_peptaibol.specificity_freq_matrix,1)
    [~,substrate_mode(i)]=max(best_peptaibol.specificity_freq_matrix(i,:));
end
best_peptaibol.substrate_mode=substrate_mode;
best_peptaibol.unique_substrate=Known_peptaibol_struct.unique_substrate;
%% 
best_peptaibol.peptaibol_list=[];
best_peptaibol.module_list=[];
best_peptaibol.seq_list=[];
for i = 1:best_peptaibol.num
    best_peptaibol.peptaibol_list=[best_peptaibol.peptaibol_list;repmat(i,best_peptaibol.A_num(i),1)];
    best_peptaibol.module_list=[best_peptaibol.module_list;[1:best_peptaibol.A_num(i)]'];
    best_peptaibol.seq_list=[best_peptaibol.seq_list;best_peptaibol.seq_matrix{i}];
end
%% 
motif_sepcial_align={'GEIVQGPTLLREYL','GEI-VQGPTLLREYL'};
ref_index=find(best_peptaibol.peptaibol_list==find(ismember(best_peptaibol.peptaibol_name,'LogPS1')));
best_peptaibol = peptaibol_dist_calculate(best_peptaibol,ref_index,motif_sepcial_align);
%% 
PS_name_raw=readcell('PS names.xlsx','Sheet','Sheet2');
%% 
for i = 1:best_peptaibol.num
    best_peptaibol.peptaibol_name{i}=PS_name_raw{ismember(PS_name_raw(:,1),best_peptaibol.peptaibol_name{i}),2};
end
%% c
best_peptaibol.peptaibol_name_str=cell(best_peptaibol.num,1);
for i = 1:best_peptaibol.num
    loc_str=split(best_peptaibol.peptaibol_name{i},'_');
    if length(loc_str)==1
        best_peptaibol.peptaibol_name_str{i}=loc_str{1};
    else
        best_peptaibol.peptaibol_name_str{i}=[loc_str{1},'_{',loc_str{2},'}'];
    end
end
%% 
ancestor=readmatrix('Peptaibol_evolution.xlsx','Sheet','Ancestor','Range','A1:U20');
best_peptaibol.module_ancestor_list=zeros(length(best_peptaibol.module_list),1);
for i = 1:best_peptaibol.num
    best_peptaibol.module_ancestor_list(best_peptaibol.peptaibol_list==i)=ancestor(ancestor(:,1)==i,2:best_peptaibol.A_num(i)+1);
end
%%
best_peptaibol.Ancestor_substrate_index=cell(max(best_peptaibol.module_ancestor_list),1);
for i =1:length(best_peptaibol.module_ancestor_list)
    if best_peptaibol.peptaibol_list(i)~=Tex2_index
        best_peptaibol.Ancestor_substrate_index{best_peptaibol.module_ancestor_list(i)}=[best_peptaibol.Ancestor_substrate_index{best_peptaibol.module_ancestor_list(i)};i];
    end
end
%% 
best_peptaibol.Ancestor_substrate_freq=[];
for i = 1:length(best_peptaibol.Ancestor_substrate_index)
    best_peptaibol.Ancestor_substrate_freq=[best_peptaibol.Ancestor_substrate_freq;mean(best_peptaibol.specificity_freq_matrix(best_peptaibol.Ancestor_substrate_index{i},:))];
end
%% 
best_peptaibol.Ancestor_substrate_entropy=zeros(length(best_peptaibol.Ancestor_substrate_index),1);
for i = 1:length(best_peptaibol.Ancestor_substrate_entropy)
    [best_peptaibol.Ancestor_substrate_entropy(i),~]=Entropy(best_peptaibol.Ancestor_substrate_freq(i,:));
end
%% 
Ancestor_module_list=[1,7,13,14,16];
for i = 1:length(Ancestor_module_list)
    loc_index=find(best_peptaibol.module_ancestor_list==Ancestor_module_list(i)&best_peptaibol.peptaibol_list~=Tex2_index);
    loc_peptaibol_list=best_peptaibol.peptaibol_list(loc_index);
    loc_module_list=best_peptaibol.module_list(loc_index);
    test_header=[];
    loc_seq=[];
    for j = 1:length(loc_peptaibol_list)
        test_header=[test_header;{[num2str(loc_index(j)),'|',num2str(loc_peptaibol_list(j)),'|',num2str(loc_module_list(j))]}];
        loc_seq=[loc_seq;best_peptaibol.A_seqs_list{loc_peptaibol_list(j)}(loc_module_list(j))];
    end
    fastawrite(['./reference_tree/peptaibol19_A_Ancestor_module_',num2str(Ancestor_module_list(i)),'.fasta'],test_header,loc_seq)
    command = ['clustalo -i ./reference_tree/peptaibol19_A_Ancestor_module_',num2str(Ancestor_module_list(i)),'.fasta -o ./reference_tree/peptaibol19_A_Ancestor_module_',num2str(Ancestor_module_list(i)),'_MSA.fasta --outfmt fasta --output-order input-order'];
    system(command);
end
%% for consensus tree
Ancestor_module_list=[1,7,13,14,16];
fileID_w = fopen('./reference_tree/peptaibol19_A_Ancestor_module_merge.treefile', 'w');
for i = 1:length(Ancestor_module_list)
    fileID = fopen(['./reference_tree/peptaibol19_A_Ancestor_module_',num2str(Ancestor_module_list(i)),'.treefile'], 'r');
    treeContent = fscanf(fileID, '%c');
    loc_seq=fastaread(['./reference_tree/peptaibol19_A_Ancestor_module_',num2str(Ancestor_module_list(i)),'.fasta']);
    for j = length(loc_seq):-1:1
        loc_header=loc_seq(j).Header;
        loc_header=split(loc_header,'|');
        if length(strfind(treeContent,loc_seq(j).Header))==1
            treeContent=strrep(treeContent,loc_seq(j).Header,loc_header{2});
        else
            j
        end
    end
    fprintf(fileID_w, '%s\n', treeContent);
end
fclose(fileID_w);
%% 
loc_seq=[];
for i = 1:best_peptaibol.num
    loc_seq=[loc_seq;best_peptaibol.A_seqs_list{i}];
end
test_header = cellstr(num2str([1:length(loc_seq)]'));
fastawrite('peptaibol19_A.fasta',test_header,loc_seq)
%% 
AA34_raw=readcell('peptaibol19_A.sig.csv');
AA34_raw(1,:)=[];
best_peptaibol.AA34=AA34_raw(:,2);
best_peptaibol.AA34_p_distance_dist=seqpdist(best_peptaibol.AA34,'ScoringMatrix','BLOSUM62','Method','p-distance','SquareForm',true);
best_peptaibol.AA34_align_score_dist=seqpdist(best_peptaibol.AA34,'ScoringMatrix','BLOSUM62','Method','alignment-score','SquareForm',true);
%% 
peptaibol_index=Module_num>=17&Contain_one_AT==1&End_with_TD_NAD>0&Module_iscomplete==1; % BGC_isfragmented==0 is selected firstly
%% 
long_peptaibol=[];% Trichoderma atrobrunneum doesn't have good BGC. All other species have at least one BGC to use.
% Trichoderma pleuroti TPhu1 need be added finally.
long_peptaibol.index=find(peptaibol_index);
long_peptaibol.Species_short=Species_short(peptaibol_index);
long_peptaibol.num=sum(peptaibol_index);
long_peptaibol.Module_num=Module_num(peptaibol_index);
long_peptaibol.peptaibol_str=cell(sum(peptaibol_index),1);
for i = 1:long_peptaibol.num
    long_peptaibol.peptaibol_str{i}=[long_peptaibol.Species_short{i},'_',num2str(long_peptaibol.Module_num(i))];
end
% sp. isn't seen as one species
[long_peptaibol.uni_peptaibol_str,~,ic]=unique(long_peptaibol.peptaibol_str);% only Trichoderma harzianum have more than one different length peptaibols
long_peptaibol.uni_species_index=cell(length(long_peptaibol.uni_peptaibol_str),1);
for i = 1:length(long_peptaibol.uni_peptaibol_str)
    long_peptaibol.uni_species_index{i}=long_peptaibol.index(ic==i);
end
loc_index=find(startsWith(long_peptaibol.uni_peptaibol_str,'Trichoderma sp.'));
for i = 1:length(loc_index)
    for j = 1:length(long_peptaibol.uni_species_index{loc_index(i)})
        long_peptaibol.uni_species_index=[long_peptaibol.uni_species_index;{long_peptaibol.uni_species_index{loc_index(i)}(j)}];
        long_peptaibol.uni_peptaibol_str=[long_peptaibol.uni_peptaibol_str;{strrep(long_peptaibol.uni_peptaibol_str{loc_index(i)},'Trichoderma sp.',Strain{long_peptaibol.uni_species_index{loc_index(i)}(j)})}];
    end
end
long_peptaibol.uni_peptaibol_str(loc_index)=[];
long_peptaibol.uni_species_index(loc_index)=[];
%% 
known_species_best_index=best_peptaibol.species_best_index;
loc_index=1:length(known_species_best_index);
loc_index(best_peptaibol.A_num==14)=[];
loc_index(known_species_best_index>length(Known_peptaibol_struct.peptaibol_index))=[];
known_species_best_index(best_peptaibol.A_num==14)=[]; % remove short peptaibol
known_species_best_index(known_species_best_index>length(Known_peptaibol_struct.peptaibol_index))=[];%remove Trichoderma pleuroti TPhu1
known_species_best_index=Known_peptaibol_struct.peptaibol_index(known_species_best_index);
long_peptaibol.known_product=zeros(length(long_peptaibol.uni_peptaibol_str),1);
long_peptaibol.known_index=zeros(length(long_peptaibol.uni_peptaibol_str),1);
long_peptaibol.species_best_index=zeros(length(long_peptaibol.uni_peptaibol_str),1);
for i = 1:length(long_peptaibol.uni_species_index)
    if any(ismember(known_species_best_index,long_peptaibol.uni_species_index{i}))
        long_peptaibol.known_product(i)=1;
        long_peptaibol.known_index(i)=loc_index(ismember(known_species_best_index,long_peptaibol.uni_species_index{i}));
        long_peptaibol.species_best_index(i)=known_species_best_index(ismember(known_species_best_index,long_peptaibol.uni_species_index{i}));
    end
end
%% 
long_peptaibol.result_list=cell(length(long_peptaibol.uni_species_index),1);
long_peptaibol.original_ntseq_list=cell(length(long_peptaibol.uni_species_index),1);
for i = 1:length(long_peptaibol.uni_species_index)
    long_peptaibol.A_num(i,1)=str2double(long_peptaibol.uni_peptaibol_str{i}(end-1:end));
    if long_peptaibol.known_index(i)==0
        long_peptaibol.original_ntseq_list{i}=cell(length(long_peptaibol.uni_species_index{i}),1);
        for j = 1:length(long_peptaibol.uni_species_index{i})
            if ~isempty(strfind(uni_Aib_locustag{long_peptaibol.uni_species_index{i}(j)},', '))% BGC_isfragmented==1. Selected firstly BGC_isfragmented==0
                loc_seq=[];
                loc_uni_Aib_locustag=split(uni_Aib_locustag{long_peptaibol.uni_species_index{i}(j)},', ');
                for k = 1:length(loc_uni_Aib_locustag)
                    loc_loc_seq=join(my_omains.seq_ntaa(ismember(my_omains.locustag,loc_uni_Aib_locustag{k})&(my_omains.region_ids==uni_region_ids(long_peptaibol.uni_species_index{i}(j))),1),'');
                    loc_seq=[loc_seq,loc_loc_seq{1}];
                end
                long_peptaibol.original_ntseq_list{i}{j}=loc_seq;
            elseif strcmp(uni_Aib_locustag{long_peptaibol.uni_species_index{i}(j)},'FUN_006869')
                loc_seq=fastaread('../data/FUN_006869.fasta');
                long_peptaibol.original_ntseq_list{i}{j}=loc_seq.Sequence;
            else
                long_peptaibol.original_ntseq_list{i}(j)=join(my_omains.seq_ntaa(ismember(my_omains.locustag,uni_Aib_locustag{long_peptaibol.uni_species_index{i}(j)})&(my_omains.region_ids==uni_region_ids(long_peptaibol.uni_species_index{i}(j))),1),'');
            end
            long_peptaibol.NT_seqs_list{i,1}{j,1}=long_peptaibol.original_ntseq_list{i}{j};
            long_peptaibol.AA_seqs_list{i,1}{j,1}=nt2aa(long_peptaibol.original_ntseq_list{i}{j},'AlternativeStartCodons',false,'ACGTOnly',false);
            X_AA_index=strfind(long_peptaibol.AA_seqs_list{i}{j},'X');
            if ~isempty(X_AA_index)
                X_NT_index=[3*X_AA_index-2,3*X_AA_index-1,3*X_AA_index];
                long_peptaibol.AA_seqs_list{i}{j}(X_AA_index)=[];
                long_peptaibol.NT_seqs_list{i}{j}(X_NT_index)=[];
            end
            loc_seq=[];
            loc_seq.Header=num2str(i);
            loc_seq.Sequence=long_peptaibol.AA_seqs_list{i}{j};
            long_peptaibol.result_list{i}{j,1}=Find_NRPS_motif_module_pfam_HRL(loc_seq,[],[],{'Aalpha','G','Talpha'},0,[],[]);
        end
    end
end
for i = 1:length(long_peptaibol.uni_species_index)
    if long_peptaibol.known_index(i)==0
        for j = 1:length(long_peptaibol.uni_species_index{i})
            assert(long_peptaibol.A_num(i)==sum(long_peptaibol.result_list{i}{j}.domain_list==2))
        end
    end
end
%% 
long_peptaibol.module_list=cell(length(long_peptaibol.uni_species_index),1);
long_peptaibol.seq_matrix=cell(length(long_peptaibol.uni_species_index),1);
long_peptaibol.NT_seq_matrix=cell(length(long_peptaibol.uni_species_index),1);
long_peptaibol.bad_motif_len_list=cell(length(long_peptaibol.uni_species_index),1);
long_peptaibol.motifid_mat=Known_peptaibol_struct.motifid_mat;
for i = 1:length(long_peptaibol.uni_species_index)
    if long_peptaibol.known_index(i)==0
        long_peptaibol.module_list{i}=[1:long_peptaibol.A_num(i)]';
        long_peptaibol.seq_matrix{i}=cell(length(long_peptaibol.uni_species_index{i}),1);
        long_peptaibol.NT_seq_matrix{i}=cell(length(long_peptaibol.uni_species_index{i}),1);
        for j = 1:length(long_peptaibol.uni_species_index{i})
            stop_index=strfind(long_peptaibol.AA_seqs_list{i}{j},'*');
            if length(stop_index)>1
                for kk = 1:length(stop_index)-1
                    long_peptaibol.AA_seqs_list{i}{j}(stop_index(kk))=[];
                    long_peptaibol.NT_seqs_list{i}{j}(3*stop_index(kk)-2:3*stop_index(kk))=[];
                end
            end
            loc_index=find(ismember(long_peptaibol.result_list{i}{j}.motifid_mat,[2,1],'rows'));
            for kk = 1:length(loc_index)
                long_peptaibol.seq_matrix{i}{j}=[long_peptaibol.seq_matrix{i}{j};[long_peptaibol.result_list{i}{j}.seq_list(loc_index(kk):loc_index(kk)+23)]'];
                loc_A_seq=join(long_peptaibol.seq_matrix{i}{j}(kk,:),'');
                loc_A_index=strfind(long_peptaibol.AA_seqs_list{i}{j},loc_A_seq{1});
                if ~isempty(loc_A_index)
                    motif_len=cellfun(@length,long_peptaibol.seq_matrix{i}{j}(kk,:));
                    for k = 1:length(motif_len)
                        long_peptaibol.NT_seq_matrix{i}{j}{kk,k}=upper(long_peptaibol.NT_seqs_list{i}{j}((loc_A_index+sum(motif_len(1:k-1)))*3-2:(loc_A_index+sum(motif_len(1:k))-1)*3));
                        loc_NT2AA=nt2aa(long_peptaibol.NT_seq_matrix{i}{j}{kk,k},'AlternativeStartCodons',false,'ACGTOnly',false);
                        if isempty(long_peptaibol.seq_matrix{i}{j}{kk,k})&&isempty(loc_NT2AA)
                            loc_NT2AA=long_peptaibol.seq_matrix{i}{j}{kk,k};
                        end
                        if ~strcmp(long_peptaibol.seq_matrix{i}{j}{kk,k},loc_NT2AA)
                            fprintf('The AA is %s, the NT is %s, the NT2AA is %s\n',long_peptaibol.seq_matrix{i}{j}{kk,k},long_peptaibol.NT_seq_matrix{i}{j}{kk,k},loc_NT2AA)
                        end
                    end
                else
                    fprintf('i=%d,j=%d,kk=%d\n',i,j,kk)
                end
            end
        end
        if length(long_peptaibol.uni_species_index{i})>1
            long_peptaibol.bad_motif_len_list{i}=zeros(length(long_peptaibol.uni_species_index{i}),1);
            for j = 1:length(long_peptaibol.uni_species_index{i})
                loc_seq_len=cellfun(@length,long_peptaibol.seq_matrix{i}{j});
                for k = 1:size(long_peptaibol.motifid_mat,1)
                    if k~=13 % except G-motif
                        if (mod(k,2)~=0&&~all(loc_seq_len(:,k)==seq_len_mode(k)))||(mod(k,2)==0&&any(loc_seq_len(:,k)==0)) %motif = mode len & intermotif ~=0
                            long_peptaibol.bad_motif_len_list{i}(j)=1;
                        end
                    end
                end
            end
        end
    end
end
%% 
for i = 1:length(Known_peptaibol_struct.NT_seqs)
    loc_index=find(Known_peptaibol_struct.species_list==i);
    for j = 1:length(loc_index)
        loc_A_seq=join(Known_peptaibol_struct.seq_matrix(loc_index(j),:),'');
        loc_A_index=strfind(Known_peptaibol_struct.AA_seqs{i},loc_A_seq{1});
        if ~isempty(loc_A_index)
            motif_len=cellfun(@length,Known_peptaibol_struct.seq_matrix(loc_index(j),:));
            for k = 1:length(motif_len)
                NT_seq_matrix{loc_index(j),k}=upper(Known_peptaibol_struct.NT_seqs{i}((loc_A_index+sum(motif_len(1:k-1)))*3-2:(loc_A_index+sum(motif_len(1:k))-1)*3));
                loc_NT2AA=nt2aa(NT_seq_matrix{loc_index(j),k},'AlternativeStartCodons',false,'ACGTOnly',false);
                if isempty(Known_peptaibol_struct.seq_matrix{loc_index(j),k})&&isempty(loc_NT2AA)
                    loc_NT2AA=Known_peptaibol_struct.seq_matrix{loc_index(j),k};
                end
                if ~strcmp(Known_peptaibol_struct.seq_matrix{loc_index(j),k},loc_NT2AA)
                    fprintf('The AA is %s, the NT is %s, the NT2AA is %s\n',Known_peptaibol_struct.seq_matrix{loc_index(j),k},NT_seq_matrix{loc_index(j),k},loc_NT2AA)
                end
            end
        else
            j % the 17th module of 73 peptaibol has problem,but 73 is dead. Known_peptaibol_struct.dead_peptaibol(73)==1. omit this problem
        end
    end
end
%% 
loc_pwd=pwd;
if exist('tmp_24601','dir')==7
    rmdir('tmp_24601', 's')
end
mkdir tmp_24601 % bulid a temporary folder
long_peptaibol.msa_matrix=cell(length(long_peptaibol.uni_species_index),1);
long_peptaibol.p_dist_matrix=cell(length(long_peptaibol.uni_species_index),1);
for i = 1:length(long_peptaibol.uni_species_index)
    if long_peptaibol.known_index(i)==0
        loc_uni_species_index=long_peptaibol.uni_species_index{i};
        loc_uni_species_index(long_peptaibol.bad_motif_len_list{i}==1)=[];
        if length(loc_uni_species_index)>1
            long_peptaibol.msa_matrix{i}=cell(long_peptaibol.A_num(i),1);
            long_peptaibol.p_dist_matrix{i}=cell(long_peptaibol.A_num(i),1);
            for j = 1:long_peptaibol.A_num(i)
                long_peptaibol.msa_matrix{i}{j}=cell(1,size(long_peptaibol.motifid_mat,1));
                for k = 1:size(long_peptaibol.motifid_mat,1)
                    if mod(k,2)~=0 %motif
                        long_peptaibol.msa_matrix{i}{j}{k}=[];
                        for kk = 1:length(long_peptaibol.uni_species_index{i})
                            long_peptaibol.msa_matrix{i}{j}{k}=[long_peptaibol.msa_matrix{i}{j}{k};long_peptaibol.seq_matrix{i}{kk}{j,k}];
                        end
                    else
                        loc_seq=[];
                        for kk = 1:length(long_peptaibol.uni_species_index{i})
                            loc_seq=[loc_seq;long_peptaibol.seq_matrix{i}{kk}(j,k)];
                        end
                        test_header=cellstr(num2str([1:length(long_peptaibol.uni_species_index{i})]'));
                        fastawrite([pwd,'/tmp_24601/',num2str(i),'_',num2str(j),'_',num2str(k),'_',num2str(kk),'.fasta'],test_header,loc_seq);
                        command = ['clustalo -i ',loc_pwd,'/tmp_24601/',num2str(i),'_',num2str(j),'_',num2str(k),'_',num2str(kk),'.fasta -o ',loc_pwd,'/tmp_24601/',num2str(i),'_',num2str(j),'_',num2str(k),'_',num2str(kk),'_1.fasta --outfmt fasta --output-order input-order'];
                        system(command);
                        raw = fastaread([pwd,'/tmp_24601/',num2str(i),'_',num2str(j),'_',num2str(k),'_',num2str(kk),'_1.fasta']);
                        long_peptaibol.msa_matrix{i}{j}{k}=[];
                        for kk = 1:length(long_peptaibol.uni_species_index{i})
                            long_peptaibol.msa_matrix{i}{j}{k}=[long_peptaibol.msa_matrix{i}{j}{k};raw(kk).Sequence];
                        end
                    end
                end
                loc_seq=[];
                for k = 1:size(long_peptaibol.motifid_mat,1)
                    loc_seq=[loc_seq,long_peptaibol.msa_matrix{i}{j}{k}];
                end
                long_peptaibol.p_dist_matrix{i}{j}=seqpdist(loc_seq,'ScoringMatrix','BLOSUM62','Method','p-distance','SquareForm',true);
            end
        end
    end
end
rmdir('tmp_24601', 's')
%% 
best_long_peptaibol=[];
best_long_peptaibol.num=length(long_peptaibol.uni_peptaibol_str);
best_long_peptaibol.uni_peptaibol_str=long_peptaibol.uni_peptaibol_str;
best_long_peptaibol.A_num=long_peptaibol.A_num;
best_long_peptaibol.species_best_index=long_peptaibol.species_best_index;
best_long_peptaibol.motifid_mat=long_peptaibol.motifid_mat;
best_long_peptaibol.seq_matrix=cell(best_long_peptaibol.num,1);%A3 A4 A5 A6
best_long_peptaibol.NT_seq_matrix=cell(best_long_peptaibol.num,1);
best_long_peptaibol.A_seqs_list=cell(best_long_peptaibol.num,1);
for i = 1:best_long_peptaibol.num
    if long_peptaibol.known_index(i)==0
        if length(long_peptaibol.uni_species_index{i})>1
            loc_uni_species_index=long_peptaibol.uni_species_index{i};
            loc_uni_species_index(long_peptaibol.bad_motif_len_list{i}==1)=[];
            if length(loc_uni_species_index)==1
                best_long_peptaibol.species_best_index(i)=loc_uni_species_index;
                best_long_peptaibol.result_list{i}=long_peptaibol.result_list{i}{ismember(long_peptaibol.uni_species_index{i},loc_uni_species_index)};
                best_long_peptaibol.NT_seqs_list{i,1}=lower(long_peptaibol.NT_seqs_list{i}{ismember(long_peptaibol.uni_species_index{i},loc_uni_species_index)});
                best_long_peptaibol.AA_seqs_list{i,1}=long_peptaibol.AA_seqs_list{i}{ismember(long_peptaibol.uni_species_index{i},loc_uni_species_index)};
                loc_seq_matrix=long_peptaibol.seq_matrix{i}{ismember(long_peptaibol.uni_species_index{i},loc_uni_species_index)};
                loc_NT_seq_matrix=long_peptaibol.NT_seq_matrix{i}{ismember(long_peptaibol.uni_species_index{i},loc_uni_species_index)};
            else
                loc_dist=sum(long_peptaibol.p_dist_matrix{i}{1});
                for j = 2:best_long_peptaibol.A_num(i)
                    loc_dist=loc_dist+sum(long_peptaibol.p_dist_matrix{i}{j});
                end
                [min_dist,I]=min(loc_dist);
                best_long_peptaibol.species_best_index(i)=long_peptaibol.uni_species_index{i}(I);
                best_long_peptaibol.result_list{i}=long_peptaibol.result_list{i}{I};
                best_long_peptaibol.NT_seqs_list{i,1}=lower(long_peptaibol.NT_seqs_list{i}{I});
                best_long_peptaibol.AA_seqs_list{i,1}=long_peptaibol.AA_seqs_list{i}{I};
                loc_seq_matrix=long_peptaibol.seq_matrix{i}{I};
                loc_NT_seq_matrix=long_peptaibol.NT_seq_matrix{i}{I};
            end
        else
            best_long_peptaibol.species_best_index(i)=long_peptaibol.uni_species_index{i};
            best_long_peptaibol.result_list{i}=long_peptaibol.result_list{i}{1};
            best_long_peptaibol.NT_seqs_list{i,1}=lower(long_peptaibol.NT_seqs_list{i}{1});
            best_long_peptaibol.AA_seqs_list{i,1}=long_peptaibol.AA_seqs_list{i}{1};
            loc_seq_matrix=long_peptaibol.seq_matrix{i}{1};
            loc_NT_seq_matrix=long_peptaibol.NT_seq_matrix{i}{1};
        end
        best_long_peptaibol.A_seqs_list{i}=join(loc_seq_matrix,'');
        for j = 1:best_long_peptaibol.A_num(i)
            for k = 7:11 %A3-A6
                best_long_peptaibol.seq_matrix{i,1}{j,k-6}=loc_seq_matrix{j,k};
                best_long_peptaibol.NT_seq_matrix{i,1}{j,k-6}=loc_NT_seq_matrix{j,k};
            end
            %merge G motif due to some sequences don't have G-motif
            loc_seq=vertcat(loc_seq_matrix(j,12:14));
            merge_seq=cell(size(loc_seq,1),1);
            for k = 1:size(loc_seq,1)
                merge_seq{k}=[loc_seq{k,1} loc_seq{k,2} loc_seq{k,3}];
            end
            best_long_peptaibol.seq_matrix{i,1}(j,6)=merge_seq;
            loc_seq=vertcat(loc_NT_seq_matrix(j,12:14));
            merge_seq=cell(size(loc_seq,1),1);
            for k = 1:size(loc_seq,1)
                merge_seq{k}=[loc_seq{k,1} loc_seq{k,2} loc_seq{k,3}];
            end
            best_long_peptaibol.NT_seq_matrix{i,1}(j,6)=merge_seq;
            for k = 15 %A3-A6
                best_long_peptaibol.seq_matrix{i,1}{j,k-8}=loc_seq_matrix{j,k};
                best_long_peptaibol.NT_seq_matrix{i,1}{j,k-8}=loc_NT_seq_matrix{j,k};
            end
        end
    else
        loc_index=find(Known_peptaibol_struct.peptaibol_index==best_long_peptaibol.species_best_index(i));
        best_long_peptaibol.result_list{i,1}=Known_peptaibol_struct.result{loc_index};
        best_long_peptaibol.NT_seqs_list{i,1}=lower(Known_peptaibol_struct.NT_seqs{loc_index});
        best_long_peptaibol.AA_seqs_list{i,1}=Known_peptaibol_struct.AA_seqs{loc_index};
        loc_index=best_peptaibol.species_best_index==loc_index;
        best_long_peptaibol.seq_matrix{i}=best_peptaibol.seq_matrix{loc_index};
        best_long_peptaibol.NT_seq_matrix{i}=best_peptaibol.NT_seq_matrix{loc_index};
        best_long_peptaibol.A_seqs_list{i}=best_peptaibol.A_seqs_list{loc_index};
    end
end
%% 
for i = 1:best_long_peptaibol.num
    if long_peptaibol.known_index(i)==0
        if length(long_peptaibol.uni_species_index{i})>1
            loc_uni_species_index=long_peptaibol.uni_species_index{i};
            loc_uni_species_index(long_peptaibol.bad_motif_len_list{i}==1)=[];
            if length(loc_uni_species_index)==1
                best_long_peptaibol.species_best_index(i)=loc_uni_species_index;
                best_long_peptaibol.result_list{i}=long_peptaibol.result_list{i}{ismember(long_peptaibol.uni_species_index{i},loc_uni_species_index)};
                best_long_peptaibol.NT_seqs_list{i,1}=long_peptaibol.NT_seqs_list{i}{ismember(long_peptaibol.uni_species_index{i},loc_uni_species_index)};
                best_long_peptaibol.AA_seqs_list{i,1}=long_peptaibol.AA_seqs_list{i}{ismember(long_peptaibol.uni_species_index{i},loc_uni_species_index)};
            else
                loc_dist=sum(long_peptaibol.p_dist_matrix{i}{1});
                for j = 2:best_long_peptaibol.A_num(i)
                    loc_dist=loc_dist+sum(long_peptaibol.p_dist_matrix{i}{j});
                end
                [min_dist,I]=min(loc_dist);
                best_long_peptaibol.species_best_index(i)=long_peptaibol.uni_species_index{i}(I);
                best_long_peptaibol.result_list{i}=long_peptaibol.result_list{i}{I};
                best_long_peptaibol.NT_seqs_list{i,1}=long_peptaibol.NT_seqs_list{i}{I};
                best_long_peptaibol.AA_seqs_list{i,1}=long_peptaibol.AA_seqs_list{i}{I};
            end
        else
            best_long_peptaibol.species_best_index(i)=long_peptaibol.uni_species_index{i};
            best_long_peptaibol.result_list{i}=long_peptaibol.result_list{i}{1};
            best_long_peptaibol.NT_seqs_list{i,1}=long_peptaibol.NT_seqs_list{i}{1};
            best_long_peptaibol.AA_seqs_list{i,1}=long_peptaibol.AA_seqs_list{i}{1};
        end
    else
        loc_index=find(Known_peptaibol_struct.peptaibol_index==best_long_peptaibol.species_best_index(i));
        best_long_peptaibol.result_list{i,1}=Known_peptaibol_struct.result{loc_index};
        best_long_peptaibol.NT_seqs_list{i,1}=Known_peptaibol_struct.NT_seqs{loc_index};
        best_long_peptaibol.AA_seqs_list{i,1}=Known_peptaibol_struct.AA_seqs{loc_index};
    end
end
%% 
assert(all(Peptaibol_check(best_long_peptaibol.species_best_index)==1)) % check choosed is peptaibol
loc_index=ismember(best_peptaibol.peptaibol_name,'PlePS/NPS1tp');
best_long_peptaibol.num=best_long_peptaibol.num+1;
best_long_peptaibol.A_num=[best_long_peptaibol.A_num;best_peptaibol.A_num(loc_index)];
best_long_peptaibol.uni_peptaibol_str=[best_long_peptaibol.uni_peptaibol_str;{[best_peptaibol.Strain_used{loc_index},'_',num2str(best_peptaibol.A_num(loc_index))]}];
best_long_peptaibol.Strain_used=Strain(best_long_peptaibol.species_best_index);
best_long_peptaibol.Strain_used=[best_long_peptaibol.Strain_used;best_peptaibol.Strain_used(loc_index)];
best_long_peptaibol.seq_matrix=[best_long_peptaibol.seq_matrix;best_peptaibol.seq_matrix(loc_index)];
best_long_peptaibol.NT_seq_matrix=[best_long_peptaibol.NT_seq_matrix;best_peptaibol.NT_seq_matrix(loc_index)];
best_long_peptaibol.A_seqs_list=[best_long_peptaibol.A_seqs_list;best_peptaibol.A_seqs_list(loc_index)];
%% 
loc_index=ismember(Known_peptaibol_struct.peptaibol_name,'PlePS/NPS1tp');
best_long_peptaibol.result_list=[best_long_peptaibol.result_list;Known_peptaibol_struct.result(loc_index)];
best_long_peptaibol.NT_seqs_list=[best_long_peptaibol.NT_seqs_list;Known_peptaibol_struct.NT_seqs(loc_index)];
best_long_peptaibol.AA_seqs_list=[best_long_peptaibol.AA_seqs_list;Known_peptaibol_struct.AA_seqs(loc_index)];
%%
best_long_peptaibol.peptaibol_list=[];
best_long_peptaibol.module_list=[];
best_long_peptaibol.seq_list=[];
for i = 1:best_long_peptaibol.num
    best_long_peptaibol.peptaibol_list=[best_long_peptaibol.peptaibol_list;repmat(i,best_long_peptaibol.A_num(i),1)];
    best_long_peptaibol.module_list=[best_long_peptaibol.module_list;[1:best_long_peptaibol.A_num(i)]'];
    best_long_peptaibol.seq_list=[best_long_peptaibol.seq_list;best_long_peptaibol.seq_matrix{i}];
end
%% 
best_long_peptaibol.seq_len_list=cellfun(@length,best_long_peptaibol.seq_list);
best_long_peptaibol.seq_len_mode=mode(best_long_peptaibol.seq_len_list);
%% 
motif_sepcial_align={'GEIVQGPTLLREYL','GEI-VQGPTLLREYL';'DFQSKARPKK','------DFQSKARPKK';'TTLPSWV','--------TTLPSWV'};
ref_index=find(best_long_peptaibol.species_best_index==Known_peptaibol_struct.peptaibol_index(best_peptaibol.species_best_index(ismember(best_peptaibol.peptaibol_name,'LogPS1'))));
best_long_peptaibol = peptaibol_dist_calculate(best_long_peptaibol,ref_index,motif_sepcial_align);
%% 
ancestor=readmatrix('Peptaibol_evolution.xlsx','Sheet','Ancestor_all','Range','A1:U40');
best_long_peptaibol.module_ancestor_list=zeros(length(best_long_peptaibol.module_list),1);
for i = 1:best_long_peptaibol.num
    best_long_peptaibol.module_ancestor_list(best_long_peptaibol.peptaibol_list==i)=ancestor(ancestor(:,1)==i,2:best_long_peptaibol.A_num(i)+1);
end
Class_group=readcell('Peptaibol_evolution.xlsx','Sheet','Ancestor_all','Range','V2:V40');
Class_group_index=ancestor(2:end,1);
best_long_peptaibol.Class_group=zeros(best_long_peptaibol.num,1);
for i = 1:length(Class_group_index)
    best_long_peptaibol.Class_group(Class_group_index(i))=str2double(strrep(Class_group{i},'Group',''));
end
%% 
Ancestor_module_list=[1,7,13,14,16];
for i = 1:length(Ancestor_module_list)
    loc_index=find(best_long_peptaibol.module_ancestor_list==Ancestor_module_list(i));
    loc_peptaibol_list=best_long_peptaibol.peptaibol_list(loc_index);
    loc_module_list=best_long_peptaibol.module_list(loc_index);
    test_header=[];
    loc_seq=[];
    for j = 1:length(loc_peptaibol_list)
        test_header=[test_header;{[num2str(loc_index(j)),'|',num2str(loc_peptaibol_list(j)),'|',num2str(loc_module_list(j))]}];
        loc_seq=[loc_seq;best_long_peptaibol.A_seqs_list{loc_peptaibol_list(j)}(loc_module_list(j))];
    end
    fastawrite(['./reference_tree/peptaibol39_A_Ancestor_module_',num2str(Ancestor_module_list(i)),'.fasta'],test_header,loc_seq)
    command = ['clustalo -i ./reference_tree/peptaibol39_A_Ancestor_module_',num2str(Ancestor_module_list(i)),'.fasta -o ./reference_tree/peptaibol39_A_Ancestor_module_',num2str(Ancestor_module_list(i)),'_MSA.fasta --outfmt fasta --output-order input-order'];
    system(command);
end
%% for consensus tree
Ancestor_module_list=[1,7,13,14,16];
fileID_w = fopen('./reference_tree/peptaibol39_A_Ancestor_module_merge.treefile', 'w');
for i = 1:length(Ancestor_module_list)
    fileID = fopen(['./reference_tree/peptaibol39_A_Ancestor_module_',num2str(Ancestor_module_list(i)),'.treefile'], 'r');
    treeContent = fscanf(fileID, '%c');
    loc_seq=fastaread(['./reference_tree/peptaibol39_A_Ancestor_module_',num2str(Ancestor_module_list(i)),'.fasta']);
    for j = length(loc_seq):-1:1
        loc_header=loc_seq(j).Header;
        loc_header=split(loc_header,'|');
        if length(strfind(treeContent,loc_seq(j).Header))==1
            treeContent=strrep(treeContent,loc_seq(j).Header,loc_header{2});
        else
            j
        end
    end
    fprintf(fileID_w, '%s\n', treeContent);
end
fclose(fileID_w);
%% check distance in ancestor module
[uni_module_ancestor,~,ic]=unique(best_long_peptaibol.module_ancestor_list);
my_edges=0:0.01:0.5;
figure('Units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(uni_module_ancestor)
    loc_A4_A5_dist_matrix=best_long_peptaibol.A4_A5_dist_matrix(ic==i,ic==i);
    upper_tri_column = [];
    if sum(ic==i)>1
        for j = 2:sum(ic==i)
            upper_tri_column=[upper_tri_column,loc_A4_A5_dist_matrix(j-1,j:end)];
        end
        subplot(4,6,i)
        histogram(upper_tri_column,my_edges)
        title({['Ancestor module ',num2str(uni_module_ancestor(i))], ['(n=',num2str(sum(ic==i)),', mean=',num2str(round(mean(upper_tri_column),2)),', max=',num2str(round(max(upper_tri_column),2)),')']})
        xlabel('Distance')
        ylabel('Count')
        set(gca, 'Fontname', 'Arial');
    end
end
%% 
loc_seq=[];
for i = 1:best_long_peptaibol.num
    loc_seq=[loc_seq;best_long_peptaibol.A_seqs_list{i}];
end
test_header = cellstr(num2str([1:length(loc_seq)]'));
fastawrite('peptaibol39_A.fasta',test_header,loc_seq)
%% 
AA34_raw=readcell('peptaibol39_A.sig.csv');
AA34_raw(1,:)=[];
best_long_peptaibol.AA34=AA34_raw(:,2);
best_long_peptaibol.AA34_p_distance_dist=seqpdist(best_long_peptaibol.AA34,'ScoringMatrix','BLOSUM62','Method','p-distance','SquareForm',true);
best_long_peptaibol.AA34_align_score_dist=seqpdist(best_long_peptaibol.AA34,'ScoringMatrix','BLOSUM62','Method','alignment-score','SquareForm',true);
%% prepare dataset file for iTOL
fileID = fopen('./reference_tree/M1/iTOL_labels_M1.txt', 'w');
fprintf(fileID, '%s\n%s\n%s\n', 'LABELS','SEPARATOR COMMA','DATA');
i=1;%M1
loc_index=find(best_long_peptaibol.module_ancestor_list==Ancestor_module_list(i));
loc_peptaibol_list=best_long_peptaibol.peptaibol_list(loc_index);
loc_module_list=best_long_peptaibol.module_list(loc_index);
for j = 1:length(loc_peptaibol_list)
    locstr=best_long_peptaibol.Strain_used{loc_peptaibol_list(j)};
    locstr=strrep(locstr,'Trichoderma','T.');
    fprintf(fileID,"%s,%s\n",[num2str(loc_index(j)),'|',num2str(loc_peptaibol_list(j)),'|',num2str(loc_module_list(j))],locstr);
end
fclose(fileID);
%%
color_list={'#8aafc9','#6EB4F6','#96cac1','#2ecaaf','#afcf78','#e7d046','#98eca5','#eab375','#ea8e83'};
fileID = fopen('./reference_tree/M1/iTOL_colors_styles_M1.txt', 'w');
fprintf(fileID, '%s\n%s\n%s\n', 'TREE_COLORS','SEPARATOR SPACE','DATA');
i=1;%M1
loc_index=find(best_long_peptaibol.module_ancestor_list==Ancestor_module_list(i));
loc_peptaibol_list=best_long_peptaibol.peptaibol_list(loc_index);
loc_module_list=best_long_peptaibol.module_list(loc_index);
for j = 1:length(loc_peptaibol_list)
    fprintf(fileID,"%s label %s bold\n",[num2str(loc_index(j)),'|',num2str(loc_peptaibol_list(j)),'|',num2str(loc_module_list(j))],color_list{best_long_peptaibol.Class_group(loc_peptaibol_list(j))});
end
fclose(fileID);
%% long and short peptaibol
peptaibol_index=Peptaibol_check==1&Module_num>=7;
peptaibol_length=Module_num(peptaibol_index);
peptaibol_foldername=foldername(peptaibol_index);
peptaibol_n=zeros(folders.num,1);
peptaibol_length_distribution=zeros(folders.num,2);%L, S
peptaibol_length_str=cell(folders.num,1);
for i = 1:folders.num
    peptaibol_n(i)=sum(ismember(peptaibol_foldername,folders.foldername{i}));
    loc_peptaibol_length=peptaibol_length(ismember(peptaibol_foldername,folders.foldername{i}));
    for j = 1:length(loc_peptaibol_length)
        if loc_peptaibol_length(j)>=17
            peptaibol_length_distribution(i,1)=peptaibol_length_distribution(i,1)+1;
        else
            peptaibol_length_distribution(i,2)=peptaibol_length_distribution(i,2)+1;
        end
    end
    peptaibol_length_str{i}=[num2str(peptaibol_length_distribution(i,1)),'_',num2str(peptaibol_length_distribution(i,2))];
end
tabulate(peptaibol_length_str)
%% 
peptaibol_length2=[];
for i = 1:length(peptaibol_length_str)
    peptaibol_length2=[peptaibol_length2;[str2double(split(peptaibol_length_str{i},'_'))]'];
end
%% 
figure
histogram(peptaibol_length)
xlabel('Peptaibol length')
ylabel('Count')
title(['Peptaibol length distribution in Trichoderma genus (n=',num2str(folders.num),')'])
saveas(gcf,'./output/figure/Peptaibol_length_distribution_in_Trichoderma_genus.svg')
%% 
product_tab=cell(folders.num,1);
product_list=[];
for i = 1:folders.num
    product_tab{i}=tabulate(my_regions.uni_sorted_product_type_str(my_regions.folderid==i));
    product_list=[product_list;product_tab{i}(:,1)];
end
uni_product=unique(product_list);
%% 
load('../tools_1025/cluster_rules/cluster_rules_table.mat');
cluster_rules_table.category{ismember(cluster_rules_table.rule,'other')}='other';
cluster_rules_table.category{ismember(cluster_rules_table.rule,'NI-siderophore')}='NI-siderophore';
cluster_rules_table.category{ismember(cluster_rules_table.rule,'NRP-metallophore')}='NRP-metallophore';
%% 
uni_product_category=cell(length(uni_product),1);
for i = 1:length(uni_product)
    loc_str=split(uni_product{i},'+');
    for j = 1:length(loc_str)
        loc_str{j}=cluster_rules_table.category{ismember(cluster_rules_table.rule,loc_str{j})};
    end
    uni_product_category(i)=join(unique(loc_str),'+');
end
%% 
product_matrix=zeros(folders.num,length(uni_product));
for i = 1:folders.num
    [~,idx]=ismember(product_tab{i}(:,1),uni_product);
    for j = 1:length(idx)
        product_matrix(i,idx(j))=product_tab{i}{j,2};
    end
end
%% 
[uni_category,~,ic]=unique(uni_product_category);
category_matrix=zeros(folders.num,length(uni_category));
for i = 1:length(uni_category)
    category_matrix(:,i)=sum(product_matrix(:,ic==i),2);
end
%% 
category_table=readcell('Category.xlsx','Sheet','category_class');
category_table(1,:)=[];
[uni_class,~,ic]=unique(category_table(:,2));
class_matrix=zeros(folders.num,length(uni_class));
for i = 1:length(uni_class)
    class_matrix(:,i)=sum(category_matrix(:,ic==i),2);
end
%% 
class_n=sum(class_matrix);
[class_n,I]=sort(class_n,'descend');
class_matrix=class_matrix(:,I);
uni_class=uni_class(I);
%% prepare dataset file for iTOL
fileID = fopen('../data/BUSCO/results/iTOL_labels.txt', 'w');
fprintf(fileID, '%s\n%s\n%s\n', 'LABELS','SEPARATOR COMMA','DATA');
for i = 1:folders.num
    if strcmp(folders.foldername{i},'Trichoderma_hypoxylon')
        fprintf(fileID,"%s,%s\n",folders.foldername{i},folders.Species{i});
    else
        fprintf(fileID,"'%s',%s\n",folders.foldername{i},folders.Species{i});
    end
end
fprintf(fileID,"'%s',%s\n",'GCA_001278495.1','Escovopsis weberi');
fclose(fileID);
%% 
color_list={'#87CEFA','#00FF00','#FF0000','#0000FF'};
fileID = fopen('../data/BUSCO/results/iTOL_colorstrip.txt', 'w');
fprintf(fileID, 'DATASET_COLORSTRIP\nSEPARATOR SPACE\nDATASET_LABEL peptaibol_num\n#num=0\nCOLOR #ffffff\n#num=1\nCOLOR #87CEFA\n#num=2\nCOLOR #00FF00\n#num=3\nCOLOR #FF0000\n#num=4\nCOLOR #0000FF\nDATA\n');
for i = 1:folders.num
    if strcmp(folders.foldername{i},'Trichoderma_hypoxylon')
        fprintf(fileID,"%s %s\n",folders.foldername{i},color_list{peptaibol_n(i)});
    else
        fprintf(fileID,"'%s' %s\n",folders.foldername{i},color_list{peptaibol_n(i)});
    end
end
fprintf(fileID,"'%s' %s\n",'GCA_001278495.1','#ffffff');
fclose(fileID);
%% 
color_list={'#ffffff','#87CEFA'};
fileID = fopen('../data/BUSCO/results/iTOL_colorstrip_L.txt', 'w');
fprintf(fileID, 'DATASET_COLORSTRIP\nSEPARATOR SPACE\nDATASET_LABEL Long_peptaibol\n#num=0\nCOLOR #ffffff\n#num=1\nCOLOR #87CEFA\nDATA\n');
for i = 1:folders.num
    if strcmp(folders.foldername{i},'Trichoderma_hypoxylon')
        fprintf(fileID,"%s %s\n",folders.foldername{i},color_list{peptaibol_length2(i,1)+1});
    else
        fprintf(fileID,"'%s' %s\n",folders.foldername{i},color_list{peptaibol_length2(i,1)+1});
    end
end
fprintf(fileID,"'%s' %s\n",'GCA_001278495.1','#ffffff');
fclose(fileID);
%% 
color_list={'#ffffff','#87CEFA','#00FF00','#FF0000'};
fileID = fopen('../data/BUSCO/results/iTOL_colorstrip_S.txt', 'w');
fprintf(fileID, 'DATASET_COLORSTRIP\nSEPARATOR SPACE\nDATASET_LABEL Short_peptaibol\n#num=0\nCOLOR #ffffff\n#num=1\nCOLOR #87CEFA\n#num=2\nCOLOR #00FF00\n#num=3\nCOLOR #FF0000\nDATA\n');
for i = 1:folders.num
    if strcmp(folders.foldername{i},'Trichoderma_hypoxylon')
        fprintf(fileID,"%s %s\n",folders.foldername{i},color_list{peptaibol_length2(i,2)+1});
    else
        fprintf(fileID,"'%s' %s\n",folders.foldername{i},color_list{peptaibol_length2(i,2)+1});
    end
end
fprintf(fileID,"'%s' %s\n",'GCA_001278495.1','#ffffff');
fclose(fileID);
%% 
known_BGC={'Trichoderma virens','Trichoderma atroviride','Trichoderma aggressivum','Trichoderma hypoxylon'};
fileID = fopen('../data/BUSCO/results/iTOL_binary_known_BGC.txt', 'w');
fprintf(fileID, 'DATASET_BINARY\nSEPARATOR COMMA\nDATASET_LABEL,known_peptaibol_BGC\nCOLOR,#ff0000\nFIELD_SHAPES,6\nFIELD_LABELS,Yes\nDATA\n');
for i = 1:folders.num
    if folders.Known_peptaibol(i)==1
        if strcmp(folders.foldername{i},'Trichoderma_hypoxylon')
            fprintf(fileID,"%s,1\n",folders.foldername{i});
        elseif ismember(folders.Short_species{i},known_BGC)
            fprintf(fileID,"'%s',1\n",folders.foldername{i});
        end
    end
end
fclose(fileID);
%% 
fileID = fopen('../data/BUSCO/results/iTOL_binary_unknown_BGC.txt', 'w');
fprintf(fileID, 'DATASET_BINARY\nSEPARATOR COMMA\nDATASET_LABEL,known_peptaibol\nCOLOR,#000000\nFIELD_SHAPES,6\nFIELD_LABELS,Yes\nDATA\n');
for i = 1:folders.num
    if folders.Known_peptaibol(i)==1
        if ~ismember(folders.Short_species{i},known_BGC)
            fprintf(fileID,"'%s',1\n",folders.foldername{i});
        end
    end
end
fclose(fileID);
%% 
fileID = fopen('../data/BUSCO/results/iTOL_multibar.txt', 'w');
fprintf(fileID, 'DATASET_MULTIBAR\nSEPARATOR COMMA\nDATASET_LABEL,Secondary metabolites statistics\nCOLOR,#2e8655\nCOLOR,#eba05F\nCOLOR,#761c77\nCOLOR,#aec1d9\nCOLOR,#f3ef7e\nCOLOR,#4063ae\nCOLOR,#e86a6a\nCOLOR,#4d6448\nFIELD_COLORS,#2e8655,#eba05F,#761c77,#aec1d9,#f3ef7e,#4063ae,#e86a6a,#4d6448\nFIELD_LABELS,NRPS,PKS,terpene,NRPS-PKS hybird,RiPP,NRPS-RiPP hybird,other hybird,other\nDATA\n');
for i = 1:folders.num
    if strcmp(folders.foldername{i},'Trichoderma_hypoxylon')
        fprintf(fileID,"%s,%d,%d,%d,%d,%d,%d,%d,%d\n",folders.foldername{i},class_matrix(i,:));
    else
        fprintf(fileID,"'%s',%d,%d,%d,%d,%d,%d,%d,%d\n",folders.foldername{i},class_matrix(i,:));
    end
end
fprintf(fileID,"'%s',%d,%d,%d,%d,%d,%d,%d,%d\n",'GCA_001278495.1',0,0,0,0,0,0,0,0);
fclose(fileID);
%% 
enzyme_list={'TqaL','TqaF','TqaM'};
Tqa_blast_raw=readcell('Tqa gene BlastP.xlsx','Range','A3:D386');
species_list=[];
index_list=[];
for i = 1:size(Tqa_blast_raw,1)
    if ~all(ismissing(Tqa_blast_raw{i,1}))&&all(ismissing(Tqa_blast_raw{i,2}))
        species_list=[species_list;{'Trichoderma_hypoxylon'}];
        index_list=[index_list;i];
    elseif ~all(ismissing(Tqa_blast_raw{i,1}))&&~all(ismissing(Tqa_blast_raw{i,2}))
        Tqa_blast_raw{i,1}=strrep(folders.Species{ismember(folders.foldername,Tqa_blast_raw(i,2))},'T.','Trichoderma');
        species_list=[species_list;Tqa_blast_raw(i,2)];
        index_list=[index_list;i];
    end
end
%% 
enzyme_count=zeros(length(species_list),3);
for i = 1:length(species_list)
    if i == length(species_list)
        loc_str=Tqa_blast_raw(index_list(i):size(Tqa_blast_raw,1),4);
    else
        loc_str=Tqa_blast_raw(index_list(i):index_list(i+1)-1,4);
    end
    for j = 1:length(enzyme_list)
        enzyme_count(i,j)=sum(ismember(loc_str,enzyme_list{j}));
    end
end
%% 
color_list={'#ffffff','#87CEFA','#00FF00','#FF0000'};
for j = 1:length(enzyme_list)
    fileID = fopen(['../data/BUSCO/results/iTOL_colorstrip_',enzyme_list{j},'.txt'], 'w');
    fprintf(fileID, 'DATASET_COLORSTRIP\nSEPARATOR SPACE\nDATASET_LABEL %s\n#num=0\nCOLOR #ffffff\n#num=1\nCOLOR #87CEFA\n#num=2\nCOLOR #00FF00\n#num=3\nCOLOR #FF0000\nDATA\n',enzyme_list{j});
    for i = 1:length(species_list)
        if strcmp(species_list{i},'Trichoderma_hypoxylon')
            fprintf(fileID,"%s %s\n",species_list{i},color_list{enzyme_count(i,j)+1});
        else
            fprintf(fileID,"'%s' %s\n",species_list{i},color_list{enzyme_count(i,j)+1});
        end
    end
    fclose(fileID);
end
%% 39 peptaibol species tree
best_long_peptaibol.foldername=cell(best_long_peptaibol.num,1);
for i = 1:best_long_peptaibol.num
    loc_str=strrep(best_long_peptaibol.Strain_used{i},'Trichoderma','T.');
    if any(ismember(folders.Species,loc_str))
        loc_index=ismember(folders.Species,loc_str);
        best_long_peptaibol.foldername{i}=folders.foldername{loc_index};
    end
end
%%
all_stain_MSA=fastaread('../data/BUSCO/results/superalignment.fa');
%% 
test_Header=[];
loc_seq=[];
for i = 1:length(all_stain_MSA)
    if ismember(all_stain_MSA(i).Header,best_long_peptaibol.foldername)||strcmp(all_stain_MSA(i).Header,'GCA_001278495.1')
        test_Header=[test_Header;{all_stain_MSA(i).Header}];
        loc_seq=[loc_seq;{all_stain_MSA(i).Sequence}];
    end
end
fastawrite('./species_tree/superalignment_39.fa',test_Header,loc_seq)
%% 
fileID = fopen('./species_tree/iTOL_labels.txt', 'w');
fprintf(fileID, '%s\n%s\n%s\n', 'LABELS','SEPARATOR COMMA','DATA');
for i = 1:best_long_peptaibol.num
    loc_str=strrep(best_long_peptaibol.Strain_used{i},'Trichoderma','T.');
    fprintf(fileID,"%s,%s\n",best_long_peptaibol.foldername{i},loc_str);
end
fprintf(fileID,"%s,%s\n",'GCA_001278495.1','Escovopsis weberi');
fclose(fileID);
%%
color_list={'#8aafc9','#6EB4F6','#96cac1','#2ecaaf','#afcf78','#e7d046','#98eca5','#eab375','#ea8e83'};
fileID = fopen('./species_tree/iTOL_colors_styles.txt', 'w');
fprintf(fileID, '%s\n%s\n%s\n', 'TREE_COLORS','SEPARATOR SPACE','DATA');
for i = 1:best_long_peptaibol.num
    fprintf(fileID,"%s label %s bold\n",best_long_peptaibol.foldername{i},color_list{best_long_peptaibol.Class_group(i)});
end
fclose(fileID);
%% MLGO
[uni_class_group,~,ic]=unique(best_long_peptaibol.Class_group);
fileID = fopen('./MLGO/input.txt', 'w');
for i = 1:length(uni_class_group)
    fprintf(fileID,'>Group%d\n',uni_class_group(i));
    loc_str=[];
    loc_index=find(ic==i,1);
    loc_ancestor_module=best_long_peptaibol.module_ancestor_list(best_long_peptaibol.peptaibol_list==loc_index);
    for j = 1:length(loc_ancestor_module)
        loc_str=[loc_str,num2str(loc_ancestor_module(j)),' '];
    end
    fprintf(fileID,'%s$\n',loc_str);
end
fclose(fileID);
%% 
count_tab=tabulate(folders.Short_species);
[~,I]=sort(cell2mat(count_tab(:,2)),'descend');
count_tab=count_tab(I,:);
%%
Tex2_index=find(ismember(best_peptaibol.peptaibol_name,'TEX2')); % only consider long PS
Aib_Vxx_index=find((sum(best_peptaibol.specificity_freq_matrix(:,ismember(best_peptaibol.unique_substrate,{'Aib';'Vxx'})),2)==1)&(best_peptaibol.specificity_freq_matrix(:,ismember(best_peptaibol.unique_substrate,'Aib'))~=1)&(best_peptaibol.specificity_freq_matrix(:,ismember(best_peptaibol.unique_substrate,'Vxx'))~=1)&best_peptaibol.peptaibol_list~=Tex2_index);
Aib_fix_index=find(best_peptaibol.specificity_freq_matrix(:,ismember(best_peptaibol.unique_substrate,'Aib'))==1&best_peptaibol.peptaibol_list~=Tex2_index);
Vxx_fix_index=find(best_peptaibol.specificity_freq_matrix(:,ismember(best_peptaibol.unique_substrate,'Vxx'))==1&best_peptaibol.peptaibol_list~=Tex2_index);
Aib_Ala_index=find((sum(best_peptaibol.specificity_freq_matrix(:,ismember(best_peptaibol.unique_substrate,{'Ala';'Aib'})),2)==1)&(best_peptaibol.specificity_freq_matrix(:,ismember(best_peptaibol.unique_substrate,'Aib'))~=1)&(best_peptaibol.specificity_freq_matrix(:,ismember(best_peptaibol.unique_substrate,'Ala'))~=1)&best_peptaibol.peptaibol_list~=Tex2_index);
Ala_fix_index=find(best_peptaibol.specificity_freq_matrix(:,ismember(best_peptaibol.unique_substrate,'Ala'))==1&best_peptaibol.peptaibol_list~=Tex2_index);
Vxx_Lxx_index=find((sum(best_peptaibol.specificity_freq_matrix(:,ismember(best_peptaibol.unique_substrate,{'Vxx';'Lxx'})),2)==1)&(best_peptaibol.specificity_freq_matrix(:,ismember(best_peptaibol.unique_substrate,'Vxx'))~=1)&(best_peptaibol.specificity_freq_matrix(:,ismember(best_peptaibol.unique_substrate,'Lxx'))~=1)&best_peptaibol.peptaibol_list~=Tex2_index);
Lxx_fix_index=find(best_peptaibol.specificity_freq_matrix(:,ismember(best_peptaibol.unique_substrate,'Lxx'))==1&best_peptaibol.peptaibol_list~=Tex2_index);
%% umap
uni_A_color=[255,255,0;126,225,168;255,55,68;108,91,255;0,0,0;211,193,179;41,110,193]/255;
loc_index={Ala_fix_index;Aib_Ala_index;Aib_fix_index;Aib_Vxx_index;Vxx_fix_index;Vxx_Lxx_index;Lxx_fix_index};
uni_A_color_matrix=[];
promiscuity_substrate_list=[];
for i = 1:length(loc_index)
    uni_A_color_matrix=[uni_A_color_matrix;repmat(uni_A_color(i,:),length(loc_index{i}),1)];
    promiscuity_substrate_list=[promiscuity_substrate_list;repmat(i,length(loc_index{i}),1)];
end
promiscuity_all_index=[Ala_fix_index;Aib_Ala_index;Aib_fix_index;Aib_Vxx_index;Vxx_fix_index;Vxx_Lxx_index;Lxx_fix_index];

A4_A5_promiscuity_MSA=[];
for i = 3:5
    A4_A5_promiscuity_MSA=[A4_A5_promiscuity_MSA,best_peptaibol.msa_list{i}(promiscuity_all_index,:)];
end
keep_col=[];% only keep different column
for i = 1:size(A4_A5_promiscuity_MSA,2)
    if length(unique(A4_A5_promiscuity_MSA(:,i)))>1
        keep_col=[keep_col,i];
    end
end
A4_A5_promiscuity_distm= seqpdist(A4_A5_promiscuity_MSA(:,keep_col),'ScoringMatrix','BLOSUM62','Method','alignment-score','SquareForm',true);
 
[reduction2, umap, clusterIdentifiers, extras]=run_umap(A4_A5_promiscuity_distm,'randomize',true);
% figure;
hold on
scatter(reduction2(:,1), reduction2(:,2), 30, uni_A_color_matrix, 'filled');
xlabel('umap1');
ylabel('umap2');
title('Ala&Aib&Vxx&Lxx (A4-A5)');
set(gca, 'Fontname', 'Arial');
hold off
% saveas(gcf,'./output/figure/A4_A5_promiscuity_umap.svg')
%% 
cluster_index={find(reduction2(:,2)>5);find(reduction2(:,1)>0);find(reduction2(:,1)>-10&reduction2(:,1)<0);find(reduction2(:,1)>-15&reduction2(:,1)<-10&reduction2(:,2)>-10&reduction2(:,2)<0);find(reduction2(:,1)<-20);find(reduction2(:,2)<-12)};
promiscuity_substrate_cluster=cell(length(cluster_index),1);
promiscuity_module_ancestor_cluster=cell(length(cluster_index),1);
for i = 1:length(cluster_index)
    promiscuity_substrate_cluster{i}=promiscuity_substrate_list(cluster_index{i});
    promiscuity_module_ancestor_cluster{i}=best_peptaibol.module_ancestor_list(promiscuity_all_index(cluster_index{i}));
end
%% heatmap
peptaibol_str={'18peptaibol','26peptaibol'};
path_pre='output/figure/substrate_promiscuity/';
for j = 1:length(peptaibol_str)
    if j == 1
        substrate_list={'Ser';'Gly';'Ala';'Aib';'Vxx';'Lxx';'Pro';'Trp';'Phe';'Tyr';'Glu';'Gln';'Glu(OMe)'};
    elseif j == 2
        substrate_list={'Ser';'Gly';'Ala';'Aib';'Vxx';'Lxx';'Lys';'Pro';'Trp';'Phe';'Tyr';'Glu';'Gln';'Glu(OMe)'};
    end
    loc_path=[path_pre,peptaibol_str{j},'/'];
    loc_edge=readtable([loc_path,'edge.xlsx']);
    loc_node=readtable([loc_path,'node.xlsx']);
    loc_matrix=zeros(length(substrate_list));
    for i = 1:length(loc_edge.source)
        loc_index1=find(ismember(substrate_list,loc_edge.source{i}));
        loc_index2=find(ismember(substrate_list,loc_edge.target{i}));
        loc_matrix(loc_index1,loc_index2)=loc_edge.weight(i);
        loc_matrix(loc_index2,loc_index1)=loc_edge.weight(i);
    end
    loc_matrix_log2=log2(loc_matrix+1);

    width = 350;   % 图形宽度（像素），可根据需要增大以容纳条形图
    height = 290;  % 图形高度（像素）
    fig = figure;
    set(fig, 'Position', [100, 100, width * 1.5, height]);  % 增加宽度以容纳右侧条形图（乘 1.5）
    % 创建 tiledlayout：1 行 2 列
    t = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    % 左侧：热图
    nexttile(1);
    h = heatmap(loc_matrix_log2);
    h.XDisplayLabels = substrate_list;
    h.YDisplayLabels = substrate_list;
    h.Title = peptaibol_str{j};  % 可选：添加标题
    hs = struct(h);
    ylabel(hs.Colorbar, 'log2(Frequence)');
    xlabel('Substrate')
    ylabel('Substrate')
    % 右侧：横向条形图（假设显示每行总和）
    nexttile(2);
    bar_data = zeros(length(substrate_list),1);
    for i = 1:length(substrate_list)
        bar_data(i)=loc_node.Size(ismember(loc_node.Names,substrate_list{i}));
    end
    b = barh(bar_data);
    ylim([0.5, length(substrate_list) + 0.5]);  % 确保 y 轴范围匹配标签    
    ax = gca;  % 获取当前轴
    set(ax, 'YDir', 'reverse');  % 反转 y 轴方向（从上到下）
    yticklabels(ax, {});
    title('A domain count')
    saveas(gcf,['./output/figure/substrate_heatmap_',peptaibol_str{j},'.svg'])
end
