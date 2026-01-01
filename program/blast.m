%% blast for Aib synthetase
sepstr='/';
toolboxpath='/storage/disk1/HRL/Project/SM_in_genome/analysis/tools_1025/';
if exist(toolboxpath,'dir')
    addpath(toolboxpath)
else
    error('Wrong path for MATLAB additional toolbox, please check\n')
end
%% 
gbk_files=dir('../data/genome/galaxy/*.gbk');
for i = 1:length(gbk_files)
	[feature_cds,~,~,~,~] = Data_GBK2CDSMisc(fullfile(gbk_files(i).folder,gbk_files(i).name));
    table_cds=struct2table(feature_cds,'AsArray',true);
    loc_str=strrep(gbk_files(i).name,'gbk','faa');
    fastawrite(['../data/genome/galaxy/',loc_str],table_cds.locus_tag,table_cds.translation);
end
%% read result
blast_result=[];
blast_result.species_list=[];
blast_result.query=[];
blast_result.match_cds=[];
blast_result.pident=[];
blast_result.length=[];
blast_result.mismatch=[];
blast_result.gapopen=[];
blast_result.qstart=[];
blast_result.qend=[];
blast_result.sstart=[];
blast_result.send=[];
blast_result.evalue=[];
blast_result.bitscore=[];
blast_result.qlen=[];
blast_result.slen=[];
blast_result.qcovhsp=[];
loc_files=dir('../data/blastp/Fungi/*.txt');
for i = 1:length(loc_files)
    loc_result=readcell(fullfile(loc_files(i).folder,loc_files(i).name),'Delimiter','\t');
    loc_file_name=loc_files(i).name;
    loc_file_name=strrep(loc_file_name,'.txt','');
    blast_result.species_list=[blast_result.species_list;repmat({loc_file_name},size(loc_result,1),1)];
    blast_result.query=[blast_result.query;loc_result(:,1)];
    blast_result.match_cds=[blast_result.match_cds;loc_result(:,3)];
    loc_result=cell2mat(loc_result(:,9:end));
    blast_result.pident=[blast_result.pident;loc_result(:,1)];
    blast_result.length=[blast_result.length;loc_result(:,2)];
    blast_result.mismatch=[blast_result.mismatch;loc_result(:,3)];
    blast_result.gapopen=[blast_result.gapopen;loc_result(:,4)];
    blast_result.qstart=[blast_result.qstart;loc_result(:,5)];
    blast_result.qend=[blast_result.qend;loc_result(:,6)];
    blast_result.sstart=[blast_result.sstart;loc_result(:,7)];
    blast_result.send=[blast_result.send;loc_result(:,8)];
    blast_result.evalue=[blast_result.evalue;loc_result(:,9)];
    blast_result.bitscore=[blast_result.bitscore;loc_result(:,10)];
    blast_result.qlen=[blast_result.qlen;loc_result(:,11)];
    blast_result.slen=[blast_result.slen;loc_result(:,12)];
    blast_result.qcovhsp=[blast_result.qcovhsp;loc_result(:,13)];
end
%% 
blast_result_T=[];
loc_result=readcell('../data/blastp/Trichoderma_hypoxylon.txt','Delimiter','\t');
blast_result_T.species_list=repmat({'Trichoderma_hypoxylon'},size(loc_result,1),1);
blast_result_T.query=loc_result(:,1);
blast_result_T.match_cds=loc_result(:,3);
loc_result=cell2mat(loc_result(:,9:end));
blast_result_T.pident=loc_result(:,1);
blast_result_T.length=loc_result(:,2);
blast_result_T.mismatch=loc_result(:,3);
blast_result_T.gapopen=loc_result(:,4);
blast_result_T.qstart=loc_result(:,5);
blast_result_T.qend=loc_result(:,6);
blast_result_T.sstart=loc_result(:,7);
blast_result_T.send=loc_result(:,8);
blast_result_T.evalue=loc_result(:,9);
blast_result_T.bitscore=loc_result(:,10);
blast_result_T.qlen=loc_result(:,11);
blast_result_T.slen=loc_result(:,12);
blast_result_T.qcovhsp=loc_result(:,13);
%% 
loc_files=dir('../data/blastp/Trichoderma/*.txt');
for i = 1:length(loc_files)
    loc_result=readcell(fullfile(loc_files(i).folder,loc_files(i).name),'Delimiter','\t');
    loc_file_name=loc_files(i).name;
    loc_file_name=strrep(loc_file_name,'.txt','');
    blast_result_T.species_list=[blast_result_T.species_list;repmat({loc_file_name},size(loc_result,1),1)];
    blast_result_T.query=[blast_result_T.query;loc_result(:,1)];
    blast_result_T.match_cds=[blast_result_T.match_cds;loc_result(:,3)];
    loc_result=cell2mat(loc_result(:,9:end));
    blast_result_T.pident=[blast_result_T.pident;loc_result(:,1)];
    blast_result_T.length=[blast_result_T.length;loc_result(:,2)];
    blast_result_T.mismatch=[blast_result_T.mismatch;loc_result(:,3)];
    blast_result_T.gapopen=[blast_result_T.gapopen;loc_result(:,4)];
    blast_result_T.qstart=[blast_result_T.qstart;loc_result(:,5)];
    blast_result_T.qend=[blast_result_T.qend;loc_result(:,6)];
    blast_result_T.sstart=[blast_result_T.sstart;loc_result(:,7)];
    blast_result_T.send=[blast_result_T.send;loc_result(:,8)];
    blast_result_T.evalue=[blast_result_T.evalue;loc_result(:,9)];
    blast_result_T.bitscore=[blast_result_T.bitscore;loc_result(:,10)];
    blast_result_T.qlen=[blast_result_T.qlen;loc_result(:,11)];
    blast_result_T.slen=[blast_result_T.slen;loc_result(:,12)];
    blast_result_T.qcovhsp=[blast_result_T.qcovhsp;loc_result(:,13)];
end
%% 
Accession_id=readcell('/storage/disk1/HRL/Project/SM_in_genome/Fungi/Fungi_genomes_clean.xlsx','Sheet','high_quality','Range','A:A');
Accession_id(1)=[];
Genus=readcell('/storage/disk1/HRL/Project/SM_in_genome/Fungi/Fungi_genomes_clean.xlsx','Sheet','high_quality','Range','BB:BB');
Genus(1)=[];
%% 
for i = 1:length(Genus)
    if ismissing(Genus{i})
        Genus{i}=[];
    end
end
%% 
blast_result_clean=blast_result;
blast_result_clean.species_list(ismember(blast_result.species_list,Accession_id(strcmp(Genus,'Trichoderma'))))=[];
blast_result_clean.query(ismember(blast_result.species_list,Accession_id(strcmp(Genus,'Trichoderma'))))=[];
blast_result_clean.match_cds(ismember(blast_result.species_list,Accession_id(strcmp(Genus,'Trichoderma'))))=[];
blast_result_clean.pident(ismember(blast_result.species_list,Accession_id(strcmp(Genus,'Trichoderma'))))=[];
blast_result_clean.length(ismember(blast_result.species_list,Accession_id(strcmp(Genus,'Trichoderma'))))=[];
blast_result_clean.mismatch(ismember(blast_result.species_list,Accession_id(strcmp(Genus,'Trichoderma'))))=[];
blast_result_clean.gapopen(ismember(blast_result.species_list,Accession_id(strcmp(Genus,'Trichoderma'))))=[];
blast_result_clean.qstart(ismember(blast_result.species_list,Accession_id(strcmp(Genus,'Trichoderma'))))=[];
blast_result_clean.qend(ismember(blast_result.species_list,Accession_id(strcmp(Genus,'Trichoderma'))))=[];
blast_result_clean.sstart(ismember(blast_result.species_list,Accession_id(strcmp(Genus,'Trichoderma'))))=[];
blast_result_clean.send(ismember(blast_result.species_list,Accession_id(strcmp(Genus,'Trichoderma'))))=[];
blast_result_clean.evalue(ismember(blast_result.species_list,Accession_id(strcmp(Genus,'Trichoderma'))))=[];
blast_result_clean.bitscore(ismember(blast_result.species_list,Accession_id(strcmp(Genus,'Trichoderma'))))=[];
blast_result_clean.qlen(ismember(blast_result.species_list,Accession_id(strcmp(Genus,'Trichoderma'))))=[];
blast_result_clean.slen(ismember(blast_result.species_list,Accession_id(strcmp(Genus,'Trichoderma'))))=[];
blast_result_clean.qcovhsp(ismember(blast_result.species_list,Accession_id(strcmp(Genus,'Trichoderma'))))=[];
%%
blast_result_all=blast_result_clean;
blast_result_all.species_list=[blast_result_all.species_list;blast_result_T.species_list];
blast_result_all.query=[blast_result_all.query;blast_result_T.query];
blast_result_all.match_cds=[blast_result_all.match_cds;blast_result_T.match_cds];
blast_result_all.pident=[blast_result_all.pident;blast_result_T.pident];
blast_result_all.length=[blast_result_all.length;blast_result_T.length];
blast_result_all.mismatch=[blast_result_all.mismatch;blast_result_T.mismatch];
blast_result_all.gapopen=[blast_result_all.gapopen;blast_result_T.gapopen];
blast_result_all.qstart=[blast_result_all.qstart;blast_result_T.qstart];
blast_result_all.qend=[blast_result_all.qend;blast_result_T.qend];
blast_result_all.sstart=[blast_result_all.sstart;blast_result_T.sstart];
blast_result_all.send=[blast_result_all.send;blast_result_T.send];
blast_result_all.evalue=[blast_result_all.evalue;blast_result_T.evalue];
blast_result_all.bitscore=[blast_result_all.bitscore;blast_result_T.bitscore];
blast_result_all.qlen=[blast_result_all.qlen;blast_result_T.qlen];
blast_result_all.slen=[blast_result_all.slen;blast_result_T.slen];
blast_result_all.qcovhsp=[blast_result_all.qcovhsp;blast_result_T.qcovhsp];
%% 
blast_result.num=length(blast_result.species_list);% 4053 Fungi
blast_result_T.num=length(blast_result_T.species_list);% 125 Trichoderma + 1 Trichoderma_hypoxylon (126)
blast_result_clean.num=length(blast_result_clean.species_list);% 4053 Fungi-22 Trichoderma (4031)
blast_result_all.num=length(blast_result_all.species_list);% 4031 Fungi + 126 Trichoderma
%%
% TqaL
% TqaM
% TqaF
enzyme_list={'TqaL','TqaF','TqaM'};
%% 
% species_clean_list=[];
% species_clean_uni_list=[];
for i = 1:length(enzyme_list)
    loc_species_list=blast_result_clean.species_list(strcmp(blast_result_clean.query,enzyme_list{i})&blast_result_clean.pident>=35&blast_result_clean.evalue<=10^(-5)&blast_result_clean.qcovhsp>=70);
%     species_clean_list=[species_clean_list;{loc_species_list}];
%     species_clean_uni_list=[species_clean_uni_list;{unique(loc_species_list)}];
    fprintf('%d homology protein of %s in %d fungal genomes (without Trichoderma)\n',length(loc_species_list),enzyme_list{i},length(unique(loc_species_list)))
end
%% 
for i = 1:length(enzyme_list)
    loc_species_list=blast_result_T.species_list(strcmp(blast_result_T.query,enzyme_list{i})&blast_result_T.pident>=35&blast_result_T.evalue<=10^(-5)&blast_result_T.qcovhsp>=70);
    fprintf('%d homology protein of %s in %d Trichoderma genomes\n',length(loc_species_list),enzyme_list{i},length(unique(loc_species_list)))
end
%% 
Taxonomy_all=readcell('/storage/disk1/HRL/Project/SM_in_genome/Fungi/Fungi_genomes_clean.xlsx','Sheet','high_quality','Range','AX:BB');
Tax_type=Taxonomy_all(1,:);
Tax_type=[{'Kingdom'},Tax_type];
Taxonomy_all(1,:)=[];
for i = 1:size(Taxonomy_all,1)
    for j = 1:size(Taxonomy_all,2)
        if ~ischar(Taxonomy_all{i,j})
            Taxonomy_all{i,j}=['Unclassified',num2str(j)];
        end
    end
end
Trichoderma_tax_index=find(ismember(Taxonomy_all(:,end),'Trichoderma'));
Trichoderma_tax=Taxonomy_all(Trichoderma_tax_index(1),:);
Taxonomy_all(Trichoderma_tax_index,:)=[];
Taxonomy_all=[Taxonomy_all;repmat(Trichoderma_tax,126,1)];
Taxonomy_all=[repmat({'Fungi'},size(Taxonomy_all,1),1),Taxonomy_all];
%% 
[Taxonomy_all_sort,I]=sortrows(Taxonomy_all);
Accession_id_raw=Accession_id;
Accession_id_raw(Trichoderma_tax_index)=[];
Accession_id_raw=[Accession_id_raw;{'Trichoderma_hypoxylon'};readcell('../data/genome/fasta/Trichoderma_fasta.xlsx','Sheet','Genome Assembly Data Report','Range','A:A')];
Accession_id_raw(ismember(Accession_id_raw,'Assembly Accession'))=[];
Accession_id_raw(ismember(Accession_id_raw,'GCA_943193675.1'))=[];
Accession_id_sort=Accession_id_raw(I);
%% set colorlist
if ~exist('assign_color.mat','file')
    ColorList=[[65,140,240;252,180,65;224,64,10;5,100,146;191,191,191;26,59,105;255,227,130;18,156,221;
                202,107,75;0,92,219;243,210,136;80,99,129;241,185,168;224,131,10;120,147,190]./255;
               [127,91,93;187,128,110;197,173,143;59,71,111;104,95,126;76,103,86;112,112,124;
                72,39,24;197,119,106;160,126,88;238,208,146]./255];
    Threshold=35;
    species_clean_list=[];
    species_clean_uni_list=[];
    for i = 1:length(enzyme_list)
        loc_species_list=blast_result_all.species_list(strcmp(blast_result_all.query,enzyme_list{i})&blast_result_all.pident>=Threshold&blast_result_all.evalue<=10^(-5)&blast_result_all.qcovhsp>=70);
        species_clean_list=[species_clean_list;{loc_species_list}];
        species_clean_uni_list=[species_clean_uni_list;{unique(loc_species_list)}];
    end
    % save species with all 3 enzymes
    species_clean_uni=species_clean_uni_list{1};
    for i = 2:length(enzyme_list)
        species_clean_uni=intersect(species_clean_uni,species_clean_uni_list{i});
    end
    % 
    Taxonomy_Aib_clean=Taxonomy_all_sort(ismember(Accession_id_sort,species_clean_uni),:);
    Accession_id_clean=Accession_id_sort(ismember(Accession_id_sort,species_clean_uni),:);
    uni_Tax_type=cell(length(Tax_type),1);
    count_list=cell(length(Tax_type),1);
    for i = 1:length(Tax_type)
        [uni_Tax_type{i},~,ic]=unique(Taxonomy_Aib_clean(:,i));
        for j = 1:length(uni_Tax_type{i})
            count_list{i}(j,1)=sum(ic==j);
            count_list{i}(j,2)=length(unique(Accession_id_clean(ic==j)));
            count_list{i}(j,3)=sum(ismember(Taxonomy_Aib_clean(:,i),uni_Tax_type{i}{j}));
        end
    end
    %
    count_threshold=20;
    Large_node_n=0;
    for i = 1:length(uni_Tax_type)
        Large_node_n=Large_node_n+sum(count_list{i}(:,2)>=count_threshold);
    end
    assign_color_flag=cell(length(uni_Tax_type),1);
    for i = 1:length(assign_color_flag)
        assign_color_flag{i}=zeros(length(uni_Tax_type{i}),1);
        assign_color_flag{i}(count_list{i}(:,2)>=20)=1;
    end
    %
    assign_color_flag_all=[];
    for i = 1:length(assign_color_flag)
        assign_color_flag_all=[assign_color_flag_all;assign_color_flag{i}];
    end
    assign_color_index=[];
    for i = 1:length(assign_color_flag_all)
        if assign_color_flag_all(i)==1
            assign_color_index=[assign_color_index;sum(assign_color_flag_all(1:i))];
        else
            assign_color_index=[assign_color_index;sum(assign_color_flag_all(1:i)==0)+sum(assign_color_flag_all)];
        end
    end
    %
    ColorList_extend=[ColorList;rand(length(assign_color_flag_all)-length(ColorList),3).*.7];
    %
    assign_color_list=cell(length(uni_Tax_type),1);
    assign_color_list{1}=ColorList_extend(assign_color_index(1:length(uni_Tax_type{1})),:);
    num_sum=length(uni_Tax_type{1});
    for i = 2:length(assign_color_list)
        assign_color_list{i}=ColorList_extend(assign_color_index(1+num_sum:length(uni_Tax_type{i})+num_sum),:);
        num_sum=num_sum+length(uni_Tax_type{i});
    end
    uni_Tax_type_match=[];
    assign_color_list_match=[];
    for i = 1:length(uni_Tax_type)
        uni_Tax_type_match=[uni_Tax_type_match;uni_Tax_type{i}];
        assign_color_list_match=[assign_color_list_match;assign_color_list{i}];
    end
    save('assign_color.mat','uni_Tax_type_match','assign_color_list_match')
else
    load('assign_color.mat')
end
%% 
Threshold_list=[35,50,70];
for Threshold_index=1:length(Threshold_list)
    Threshold=Threshold_list(Threshold_index);
    species_clean_list=[];
    species_clean_uni_list=[];
    for i = 1:length(enzyme_list)
        loc_species_list=blast_result_all.species_list(strcmp(blast_result_all.query,enzyme_list{i})&blast_result_all.pident>=Threshold&blast_result_all.evalue<=10^(-5)&blast_result_all.qcovhsp>=70);
        species_clean_list=[species_clean_list;{loc_species_list}];
        species_clean_uni_list=[species_clean_uni_list;{unique(loc_species_list)}];
        fprintf('%d homology protein of %s in %d fungal genomes\n',length(loc_species_list),enzyme_list{i},length(unique(loc_species_list)))
    end
    % save species with all 3 enzymes
    species_clean_uni=species_clean_uni_list{1};
    for i = 2:length(enzyme_list)
        species_clean_uni=intersect(species_clean_uni,species_clean_uni_list{i});
    end
    % 
    Taxonomy_Aib_clean=Taxonomy_all_sort(ismember(Accession_id_sort,species_clean_uni),:);
    Accession_id_clean=Accession_id_sort(ismember(Accession_id_sort,species_clean_uni),:);
    uni_Tax_type=cell(length(Tax_type),1);
    count_list=cell(length(Tax_type),1);
    for i = 1:length(Tax_type)
        [uni_Tax_type{i},~,ic]=unique(Taxonomy_Aib_clean(:,i));
        for j = 1:length(uni_Tax_type{i})
            count_list{i}(j,1)=sum(ic==j);
            count_list{i}(j,2)=length(unique(Accession_id_clean(ic==j)));
            count_list{i}(j,3)=sum(ismember(Taxonomy_Aib_clean(:,i),uni_Tax_type{i}{j}));
        end
    end
    links=[];
    for i = 1:length(Tax_type)-1
        loc_tax=cell(length(Taxonomy_Aib_clean),1);
        for j = 1:length(Taxonomy_Aib_clean)
            loc_tax(j)=join(Taxonomy_Aib_clean(j,i:i+1),';');
        end
        [loc_tax_pair,~,ic]=unique(loc_tax,'stable');
        for j = 1:length(loc_tax_pair)
            loc_str=split(loc_tax_pair{j},';');
            loc_str{1}=[loc_str{1},'(',num2str(count_list{i}(ismember(uni_Tax_type{i},loc_str{1}),2)),')'];
            loc_str{2}=[loc_str{2},'(',num2str(count_list{i+1}(ismember(uni_Tax_type{i+1},loc_str{2}),2)),')'];
            links=[links;[loc_str',{sum(ic==j)}]];
        end
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

    width = 1300;  % Width of the figure in pixels
    height = 1000*35/Threshold;  % Height of the figure in pixels

    fig = figure;
    set(fig, 'Position', [100, 100, width, height]);  % Set the position and size

    SK=SSankey(links(:,1),links(:,2),links(:,3),'RenderingMethod','left','ColorList',loc_color_list);

    SK.draw()
    title(['Taxonomy distribution of genomes with AibABC cluster (Threshold=',num2str(Threshold),')'])
    % saveas(gcf,['./output/figure/Sankey_AibABC_distribution_T',num2str(Threshold),'.svg'])
end