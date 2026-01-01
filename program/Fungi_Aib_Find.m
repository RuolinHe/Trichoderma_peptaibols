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
    load('./output/Fungi/my_omains.mat');
    load('./output/Fungi/folders.mat');
    load('./output/simple_regions.mat')
    my_regions=simple_regions;
    fprintf('Finsh\n')
    toc
end
if ~exist('./output/figure','dir')
    mkdir('./output/figure');
end
A_domain=2;
%% 
fprintf('Find A domain motif...\t')
A_domain_index=find(my_omains.typeid_mat(:,1)==A_domain&my_omains.isdomain==1);
A_domain_folder_ids=my_regions.folderid(my_omains.region_ids(A_domain_index));
if ~exist('./output/result.mat','file')
    result=cell(length(A_domain_index),1);
    for i = 1:length(A_domain_index)
        seq=[];
        seq.Header=['omain_ids|',num2str(A_domain_index(i))];
        seq.Sequence=[my_omains.seq_ntaa{A_domain_index(i)-1,2},my_omains.seq_ntaa{A_domain_index(i),2},my_omains.seq_ntaa{A_domain_index(i)+1,2}];
        if ~isempty(seq.Sequence)
            result{i} = Find_NRPS_motif_module_pfam_HRL(seq,[],[],{'Aalpha','G','Talpha'},0,[],[]);
        end
    end
    save('./output/result','result','A_domain_index','A_domain_folder_ids','-v7.3')
else
    load('./output/result.mat')
end
fprintf('Finsh\n')
%% 
tic
fprintf('Find Aib...\t')
msa_outputpath='./output/msasummary/';
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
    save('./output/result','result','A_domain_index','A_domain_folder_ids','-v7.3')
    fprintf('Finsh\n')
end
%% 
Dist_matrix=NaN(length(result),2);
for i = 1:length(result)
    if result{i}.good==1
        Dist_matrix(i,1)=result{i}.MinDist;
        Dist_matrix(i,2)=result{i}.MeanDist;
    end
end
%% 
figure('Units','normalized','outerposition',[0 0 1 1]);
histogram(Dist_matrix(:,1))
xlabel('Min distance to Aib')
ylabel('Count')
title(['A domain min distance distribution in the ',num2str(folders.num),' fungal genomes (n=',num2str(sum(~isnan(Dist_matrix(:,1)))),')'])
saveas(gcf,'./output/figure/A_min_dist_distribution_fungi.svg')
%% 
figure('Units','normalized','outerposition',[0 0 1 1]);
histogram(log10(Dist_matrix(:,1)))
xlabel('Min distance to Aib (log10)')
ylabel('Count')
title(['A domain min distance distribution in the ',num2str(folders.num),' fungal genomes (n=',num2str(sum(~isnan(Dist_matrix(:,1)))),')'])
xlim([-2,inf])
saveas(gcf,'./output/figure/A_min_dist_distribution_fungi_log10.svg')
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
saveas(gcf,'./output/figure/A_min_dist_distribution_fungi_zoom_in.svg')
%% 
threshold=0.075;
tabulate(folders.Genus(A_domain_folder_ids(Dist_matrix(:,1)<threshold)))
%% 
[uni_region_ids,~,ic]=unique(my_omains.region_ids(A_domain_index(Dist_matrix(:,1)<threshold)));
Aib_num=zeros(length(uni_region_ids),1);% the number of A domain which specify Aib
for i = 1:length(uni_region_ids)
    Aib_num(i)=sum(ic==i);
end
uni_sorted_product_type_str=my_regions.uni_sorted_product_type_str(uni_region_ids);
%% generate a simple region with basic information
% simple_regions=[];
% simple_regions.folderid=my_regions.folderid;
% simple_regions.foldername=my_regions.foldername;
% simple_regions.filename=my_regions.filename;
% simple_regions.contig_edge=my_regions.contig_edge;
% simple_regions.product_single_type=my_regions.product_single_type;
% simple_regions.num=my_regions.num;
% simple_regions.uni_sorted_product_type_str=cell(simple_regions.num,1);
% for i = 1:simple_regions.num
%     simple_regions.uni_sorted_product_type_str(i)=join(simple_regions.product_single_type{i},'+');
% end
% save('./output/simple_regions','simple_regions')
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
Phylum=folders.Phylum(my_regions.folderid(uni_region_ids));
Class=folders.Class(my_regions.folderid(uni_region_ids));
Order=folders.Order(my_regions.folderid(uni_region_ids));
Family=folders.Family(my_regions.folderid(uni_region_ids));
Genus=folders.Genus(my_regions.folderid(uni_region_ids));
Species=folders.Species(my_regions.folderid(uni_region_ids));
Strain=folders.Strain(my_regions.folderid(uni_region_ids));
Species_short=folders.Species_short(my_regions.folderid(uni_region_ids));
foldername=my_regions.foldername(uni_region_ids);
filename=my_regions.filename(uni_region_ids);
for i = 1:length(foldername)
    loc_foldername=split(foldername{i},'/');
    foldername(i)=loc_foldername(end);
end
%%
Aib_infor=table(Species_short,Aib_num,Module_num,Start_with_Aib,AT_before_start_Aib,Contain_one_AT,End_with_PP_TD,End_with_PP_NAD,End_with_PP_TD_NAD,End_with_TD,End_with_NAD,End_with_TD_NAD,uni_sorted_product_type_str,BGC_architecture,Aib_locustag_A_ids,uni_Aib_locustag,foldername,filename,Phylum,Class,Order,Family,Genus,Species,Strain,Module_iscomplete,BGC_isfragmented);
writetable(Aib_infor,'Aib_infor.xlsx','Sheet','infor')
%% 
fprintf('Genus distribution in %d BGCs:\n',length(Genus))
tabulate(Genus)
fprintf('NRPS ending with TD domain:\n')
tabulate(Genus(End_with_PP_TD==1))
fprintf('NRPS ending with TD domain and the number of A domain with Aib substrate > 1:\n')
tabulate(Genus(End_with_PP_TD==1&Aib_num>1))
fprintf('NRPS ending with TD domain and the number of A domain > 10:\n')
tabulate(Genus(End_with_PP_TD==1&Module_num>10))
fprintf('NRPS ending with TD domain and the number of A domain with Aib substrate > 1 and the number of A domain > 10:\n')
tabulate(Genus(End_with_PP_TD==1&Aib_num>1&Module_num>10))
%% 
A_Genus=folders.Genus(A_domain_folder_ids);
A_folder=folders.foldername(A_domain_folder_ids);
clean_Dist_matrix=Dist_matrix(~ismember(A_Genus,'Trichoderma'),:);
clean_A_folder=A_folder(~ismember(A_Genus,'Trichoderma'));
uni_clean_A_folder=unique((clean_A_folder(clean_Dist_matrix(:,1)<threshold)));
sum(clean_Dist_matrix(:,1)<0.075) % number of Aib-A not in Trichoderma
%% 
load('output/result_T.mat','result')
Dist_matrix=NaN(length(result),2);
for i = 1:length(result)
    if result{i}.good==1
        Dist_matrix(i,1)=result{i}.MinDist;
        Dist_matrix(i,2)=result{i}.MeanDist;
    end
end
%% 
clean_Dist_matrix=[clean_Dist_matrix;Dist_matrix];
%% 
figure('Units','normalized','outerposition',[0 0 1 1]);
histogram(log10(clean_Dist_matrix(:,1)))
xlabel('Min distance to Aib (log10)')
ylabel('Count')
title(['A domain min distance distribution in the ',num2str(folders.num-22+126),' fungal genomes (n=',num2str(sum(~isnan(clean_Dist_matrix(:,1)))),')'])
xlim([-2,inf])
saveas(gcf,'./output/figure/A_min_dist_distribution_clean_fungi_log10.svg')
%% 
edge=0:0.005:0.15;
figure
histogram(clean_Dist_matrix(clean_Dist_matrix(:,1)<=0.15,1),edge)
hold on
xlabel('Min distance to Aib')
ylabel('Count')
loc_ylim=ylim;
line([0.075,0.075],loc_ylim,'Color','r','LineStyle','--')
hold off
saveas(gcf,'./output/figure/A_min_dist_distribution_clean_fungi_zoom_in.svg')
%% tabulate taxonomy
Taxonomy_all=readcell('../data/Fungi_genomes_clean.xlsx','Sheet','high_quality','Range','AX:BB');
Tax_type=Taxonomy_all(1,:);
Taxonomy_all(1,:)=[];
for i = 1:size(Taxonomy_all,1)
    for j = 1:size(Taxonomy_all,2)
        if ~ischar(Taxonomy_all{i,j})
            Taxonomy_all{i,j}='Unclassified';
        end
    end
end
Trichoderma_tax_index=find(ismember(Taxonomy_all(:,end),'Trichoderma'));
Trichoderma_tax=Taxonomy_all(Trichoderma_tax_index(1),:);
Taxonomy_all(Trichoderma_tax_index,:)=[];
Taxonomy_all=[Taxonomy_all;repmat(Trichoderma_tax,126,1)];
%% 
Module_num=readmatrix('Aib_infor.xlsx','Sheet','infor0','Range','E2:E379');
Peptaibol_check=readmatrix('Aib_infor.xlsx','Sheet','infor0','Range','AE2:AE379');
Accession=readcell('Aib_infor.xlsx','Sheet','infor0','Range','T2:T379');
Taxonomy_pep=readcell('Aib_infor.xlsx','Sheet','infor0','Range','V2:Z379');
Species_pep=readcell('Aib_infor.xlsx','Sheet','infor0','Range','AA2:AA379');
for i = 1:size(Taxonomy_pep,1)
    for j = 1:size(Taxonomy_pep,2)
        if ~ischar(Taxonomy_pep{i,j})
            Taxonomy_pep{i,j}='Unclassified';
        end
    end
end
Taxonomy_pep_clean=Taxonomy_pep(Peptaibol_check==1&Module_num>=7,:);
Species_pep_clean=Species_pep(Peptaibol_check==1&Module_num>=7,:);
Module_num_clean=Module_num(Peptaibol_check==1&Module_num>=7,:);
Accession_clean=Accession(Peptaibol_check==1&Module_num>=7,:);
[~,ia,~]=unique(Accession(Peptaibol_check==1&Module_num>=7,:));
Taxonomy_pep_clean_uni=Taxonomy_pep_clean(ia,:);
%% 
for i = 1:length(Tax_type)
    fprintf('%s\n',Tax_type{i});
    tabulate(Taxonomy_pep_clean(:,i))
end
%% 
for i = 1:length(Tax_type)
    fprintf('%s\n',Tax_type{i});
    tabulate(Taxonomy_pep_clean_uni(:,i))
end
%% 
uni_Tax_type=cell(length(Tax_type),1);
count_list=cell(length(Tax_type),1);
for i = 1:length(Tax_type)
    [uni_Tax_type{i},~,ic]=unique(Taxonomy_pep_clean(:,i));
    for j = 1:length(uni_Tax_type{i})
        count_list{i}(j,1)=sum(ic==j);
        count_list{i}(j,2)=length(unique(Accession_clean(ic==j)));
        count_list{i}(j,3)=sum(ismember(Taxonomy_all(:,i),uni_Tax_type{i}{j}));
    end
end
%% 
fig = figure('Units','normalized','outerposition',[0 0 1 1]);
left_color = [0 0 0];
right_color = [239/255 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
subplot(2,3,1)
yyaxis left
bar(count_list{1}(:,2),'LineStyle','none','BarWidth',0.8,'FaceColor','#737373')
xticklabels('Fungi')
ylabel('Count')
yyaxis right
scatter(1,100*count_list{1}(:,2)./4157,72,'^')
ytickformat('%.0f%%')
ylabel('Density')
ylim([0,105])
title('Kingdom')
for i = 1:length(uni_Tax_type)
    loc_x=1:length(uni_Tax_type{i});
    subplot(2,3,i+1)
    yyaxis left
    bar(loc_x,count_list{i}(:,2),'LineStyle','none','BarWidth',0.8,'FaceColor','#737373')
    xticks(loc_x)
    xticklabels(uni_Tax_type{i})
    ylabel('Count')
    yyaxis right
    scatter(loc_x,100*count_list{i}(:,2)./count_list{i}(:,3),72,'^')
    xticks(loc_x)
    ytickformat('%.0f%%')
    ylabel('Density')
    ylim([0,105])
    title(Tax_type{i})
end
saveas(gcf,'./output/figure/Bar_PS_distribution.svg')
%% 
if ~exist('Taxonomy_pep_clean_uni','var')
    load('SSankey_Fungi_Aib_Find.mat')
end
links=[];
for i = 1:length(Tax_type)-1
    loc_tax=cell(length(Taxonomy_pep_clean_uni),1);
    for j = 1:length(Taxonomy_pep_clean_uni)
        loc_tax(j)=join(Taxonomy_pep_clean_uni(j,i:i+1),';');
    end
    [loc_tax_pair,~,ic]=unique(loc_tax,'stable');
    for j = 1:length(loc_tax_pair)
        loc_str=split(loc_tax_pair{j},';');
        loc_str{1}=[loc_str{1},'(',num2str(count_list{i}(ismember(uni_Tax_type{i},loc_str{1}),2)),')'];
        loc_str{2}=[loc_str{2},'(',num2str(count_list{i+1}(ismember(uni_Tax_type{i+1},loc_str{2}),2)),')'];
        links=[links;[loc_str',{sum(ic==j)}]];
    end
end
links=[links;links(12:13,:);links(10,:)];
links([10,12:13],:)=[];

if ~exist('assign_color_list_match','var')
    load('assign_color.mat')
end
NodeList=[links(:,1);links(:,2)];
NodeList=unique(NodeList,'stable');
loc_color_list=[];
for i = 1:length(NodeList)
    loc_node=split(NodeList{i},'(');
    loc_node=loc_node{1};
    if any(contains(uni_Tax_type_match,loc_node))
        loc_color_list=[loc_color_list;assign_color_list_match(contains(uni_Tax_type_match,loc_node),:)];
    else
        error('Species is missing in uni_Tax_type_match')
    end
end
%% 
% 创建桑基图对象(Create a Sankey diagram object)
width = 1240;  % Width of the figure in pixels
height = 420;  % Height of the figure in pixels

% Create a figure with the specified size
fig = figure;
set(fig, 'Position', [100, 100, width, height]);  % Set the position and size
SK=SSankey(links(:,1),links(:,2),links(:,3),'RenderingMethod','left','ColorList',loc_color_list);

% 开始绘图(Start drawing)
SK.draw()
% saveas(gcf,'./output/figure/Sankey_PS_distribution.svg')
%% 
peptaibol_n_matrix=zeros(2);% L_strain,L_species;S_strain,S_species
loc_species=Species_pep_clean(Module_num_clean>=17);
loc_Accession_clean=Accession_clean(Module_num_clean>=17);
loc_genus=Taxonomy_pep_clean(Module_num_clean>=17,end);
[~,ia,~]=unique(loc_Accession_clean);
loc_species=loc_species(ia);
loc_genus=loc_genus(ia);
sp_n=sum(ismember(loc_species,'sp.'));
peptaibol_n_matrix(1,1)=length(loc_species);
peptaibol_n_matrix(1,2)=length(unique(loc_species));
if sp_n>0
    peptaibol_n_matrix(1,2)=peptaibol_n_matrix(1,2)+sp_n-1;
end

loc_species=Species_pep_clean(Module_num_clean<17);
loc_Accession_clean=Accession_clean(Module_num_clean<17);
[~,ia,~]=unique(loc_Accession_clean);
loc_species=loc_species(ia);
sp_n=sum(ismember(loc_species,'sp.'));
peptaibol_n_matrix(2,1)=length(loc_species);
peptaibol_n_matrix(2,2)=length(unique(loc_species));
if sp_n>0
    peptaibol_n_matrix(2,2)=peptaibol_n_matrix(1,2)+sp_n-1;
end
%% 
figure
b=bar(peptaibol_n_matrix');
ylabel('Count')
xticklabels({'Genomes','Strain'})
legend({'Long','Short'})
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom')
saveas(gcf,'./output/figure/L_S_PS_distribution.svg')