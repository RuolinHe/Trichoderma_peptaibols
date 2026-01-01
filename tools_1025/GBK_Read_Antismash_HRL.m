function rawset=GBK_Read_Antismash_HRL(datapaths,fname,sepstr,printprogress,antismashversion,regions_flag,omains_flag,CDSs_flag,cand_clusters_flag,comparison_data,merge_threshold)
% GBK_Read_Antismash_HRL is used for reading gbk files in antismash result
% The simplest usage: rawset=GBK_Read_Antismash_HRL(datapaths,fname)
% datapaths: the mian folder where the antismash result folders are stored
% fname: single antismash result folder
% Example structure of datapaths and fname
% datapaths
% └── fname
%     ├── css
%     ├── images
%     ├── svg
%     ├── knownclusterblast
%     │   └── region1
%     ├── fname.gbk
%     ├── fname.json
%     ├── xxx.region001.gbk
%     └── xxx.region002.gbk
% Note: there should be a gbk and json with the same name to fname.
% Therefore, your minimal antismash command should be
% antismash fname.gbff --output-dir ./fname
% if gname is "GCF_009498275.1", the command will be "antismash GCF_009498275.1.gbff --output-dir ./GCF_009498275.1"
% The recommended command is "antismash GCF_009498275.1.gbff --asf --clusterhmmer --cc-mibig --tfbs --cb-knownclusters --output-dir ./GCF_009498275.1 --skip-zip-file --taxon bacteria --genefinding-tool none"
% "--cb-knownclusters" is optional
%
% Other parameters:
% sepstr: it depends on your system. '\' for windows, '/' for non-windows. Default is determined by "ispc" function.
% printprogress: 0 or 1. For 1, it will print gbkfilename in processing. Default is 0.
% antismashversion: integer. Which version of antismash is used. Default is 7. And this function only test in antismash v7.0.0
% regions_flag: 0 or 1. For 1, it will generate regions struct. Default is 1.
% omains_flag: 0 or 1. For 1, it will generate omains struct. Default is 1.
% CDSs_flag: 0 or 1. For 1, it will generate CDSs struct. Default is 1.
% cand_clusters_flag: 0 or 1. For 1, it will generate cand_clusters struct. Default is 1.
% Note: product_after_merge in regions is dependent on cand_clusters. And
% omains is dependent on CDSs. It recommends that setting these flag all to
% 1 for generate completed information.
% merge_threshold: float in 0~1. If the neighboring overlap ratio exceeds this threshold, the products are merged. Default is 0.85 (test in mibig v3).
% comparison_data: the database used in cluster comprison, default is MIBiG
%%
loc_folderpath=[datapaths,fname];
if nargin<11
    merge_threshold=0.85;
end
if nargin<10
    % comparison_data='custom_db';
    comparison_data='MIBiG';
end
if nargin<5
    antismashversion=7; % this function only test in antismash v7.0
    regions_flag=1; % only save regions
    omains_flag=1;
    CDSs_flag=1;
    cand_clusters_flag=1;
end
if nargin<4
    printprogress=0;
end
if nargin<3
    if ispc
        sepstr = '\';  % Windows系统使用反斜杠
    else
        sepstr = '/';   % 非Windows系统使用正斜杠
    end
end
%% the output of the function
rawset=[];
regions=[];regions.num=0; % previously cluster or pathway, but it may not be an good cluster
CDSs=[];CDSs.num=0;
omains=[];omains.num=0;
cand_clusters=[];cand_clusters.num=0;
%%%%%%%%%%%%
filesinfolder=dir([loc_folderpath,sepstr,'*.region*.gbk']);
if ~isempty(filesinfolder)
    %% load motif information used in TFBS Finder
    % jsonData = fileread('pwms.json'); % download from https://github.com/HAugustijn/MiniMotif/blob/main/bin/data/pwms.json
    % TFBSjson = jsondecode(jsonData);
    % TFBSinfor=[];
    % TFBSinfor.name=fieldnames(TFBSjson);
    % for i = 1:length(TFBSinfor.name)
    %     TFBSinfor.description{i,1}=TFBSjson.(TFBSinfor.name{i}).description;
    %     TFBSinfor.consensus{i,1}=TFBSjson.(TFBSinfor.name{i}).consensus;
    %     TFBSinfor.max_score(i,1)=TFBSjson.(TFBSinfor.name{i}).max_score;
    %     TFBSinfor.min_score(i,1)=TFBSjson.(TFBSinfor.name{i}).min_score;
    %     TFBSinfor.pwm{i,1}=TFBSjson.(TFBSinfor.name{i}).pwm;
    % end
    %% define pattern
    resist_pattern = 'resistance \(resist\) ([\w-]+) \(Score: (\d+(\.\d+)?); E-value: ([\d.e-]+)\)';
    TFBS_pattern = 'TFBS match to (?<TFBS>\w+), (?<Description>.*?), confidence: (?<Confidence>\w+), score: (?<Score>[\d.]+)';
    region_pattern = '\.region0*(\d+)';
    t2pks_elongation_pattern='(\d+(\|\d+)*) \(Score: (\d+(\.\d+)?); E-value: ([\de.-]+)\)';
    t2pks_starter_units_pattern='([^()]+) \(Score: ([\d.]+); E-value: ([\de.-]+)\)';
    %% load domain names
    domainnames=[];
    raw=readtable('Antismash_domainnames.xlsx');% load all domain types
    domainnames.abr=raw.domainname_short;
    domainnames.formal=raw.domainname;
    domainnames.explain=raw.domainname_explain;
    domainnames.typelist=[];
    % compensate the abrivation name
    domainnames.dtype=zeros(size(domainnames.abr));
    for j=1:length(domainnames.abr)
        if strcmp(domainnames.abr{j},'')
            domainnames.abr{j}=domainnames.formal{j};
        end
        if j==1
            domainnames.dtype(j)=1;
            domainnames.typelist{1}=domainnames.abr{j};
        elseif ~strcmp(domainnames.abr{j},domainnames.abr{j-1})
            domainnames.dtype(j)=domainnames.dtype(j-1)+1;
            domainnames.typelist{domainnames.dtype(j-1)+1}=domainnames.abr{j};
        else
            domainnames.dtype(j)=domainnames.dtype(j-1);
        end
    end
    omains.domaintypelist=domainnames.typelist;

    jsonData = fileread([loc_folderpath,sepstr,fname,'.json']);
    jsonData = jsondecode(jsonData);
    % prepare for RiPP information
    ripp_str_list={'lanthipeptide';'lassopeptide';'sactipeptide';'thiopeptide'};
    ripp_list=cell(4,2);
    for i = 1:length(jsonData.records)
        if ~isempty(jsonData.records(i).areas)
            for j = 1:length(ripp_str_list)
                loc_ripp=jsonData.records(i).modules.(['antismash_modules_',ripp_str_list{j},'s']).motifs;
                if ~isempty(loc_ripp)
                    loc_fieldnames=fieldnames(loc_ripp);
                    if ~isempty(loc_fieldnames)
                        if ~strcmp(ripp_str_list{j},'thiopeptide')
                            for k = 1:length(loc_fieldnames)
                                ripp_list{j,1}=[ripp_list{j,1};repmat(loc_fieldnames(k),length(loc_ripp.(loc_fieldnames{k})),1)];%locus_tag
                                ripp_list{j,2}=[ripp_list{j,2};loc_ripp.(loc_fieldnames{k})];
                            end
                        else
                            for k = 1:length(loc_ripp)
                                ripp_list{j,1}=[ripp_list{j,1};{strrep(loc_ripp(k).locus_tag,['_',ripp_str_list{j}],'')}];%locus_tag
                                ripp_list{j,2}=[ripp_list{j,2};loc_ripp(k)];
                            end
                        end
                    end
                end
            end
        end
    end
    for j = 1:length(ripp_str_list)
        if ~isempty(ripp_list{j,1})
            for k = 1:length(ripp_list{j,2})
                if strcmp(ripp_str_list{j},'sactipeptide')
                    ripp_list{j,2}(k).RODEO_score=ripp_list{j,2}(k).score;
                    ripp_list{j,2}(k).Cleavage_pHMM_score=[];
                    ripp_list{j,2}(k).detailed_info=[];
                else
                    ripp_list{j,2}(k).RODEO_score=str2double(ripp_list{j,2}(k).detailed_info.RODEO_score{1});
                    ripp_list{j,2}(k).Cleavage_pHMM_score=ripp_list{j,2}(k).score;
                    ripp_list{j,2}(k).detailed_info=rmfield(ripp_list{j,2}(k).detailed_info,'RODEO_score');
                end
            end
            ripp_list{j,2}=rmfield(ripp_list{j,2},'score');
        end
    end
    ripp_all_list=cell(2,1);
    for i = 1:size(ripp_list,1)
        ripp_all_list{1}=[ripp_all_list{1};ripp_list{i,1}];
        ripp_all_list{2}=[ripp_all_list{2};ripp_list{i,2}];
    end
    [~,whole_feature_other,~,~,wholegbk] = Data_GBK2CDSMisc([loc_folderpath,sepstr,fname,'.gbk']);
    wholegbk_LocusName=cell(length(wholegbk),1);
    wholegbk_LocusSequenceLength=zeros(length(wholegbk),1);
    wholegbk_LocusTopology=cell(length(wholegbk),1);
    for i = 1:length(wholegbk)
        wholegbk_LocusName{i,1}=wholegbk(i).LocusName;
        wholegbk_LocusSequenceLength(i,1)=str2double(wholegbk(i).LocusSequenceLength);
        wholegbk_LocusTopology{i,1}=wholegbk(i).LocusTopology;
    end
    whole_table_other=struct2table(whole_feature_other,'AsArray',true);
    if cand_clusters_flag
        cand_clusters_id=find(ismember(whole_table_other.Category,'cand_cluster'));
        cand_clusters.num=length(cand_clusters_id);
        if isfield(whole_feature_other,'SMILES')
            cand_clusters.SMILES=whole_table_other.SMILES(cand_clusters_id);
        else
            cand_clusters.SMILES=cell(cand_clusters.num,1);
        end
        cand_clusters.border=[whole_table_other.Minimum(cand_clusters_id),whole_table_other.Maximum(cand_clusters_id)];% border in the whole genome
        cand_clusters.candidate_cluster_number=str2double(whole_table_other.candidate_cluster_number(cand_clusters_id));% in the genome fragment such as contig
        for i =1:cand_clusters.num
            cand_clusters.protoclusters{i,1}=str2double(whole_table_other.protoclusters{cand_clusters_id(i)});% in the genome fragment such as contig
        end
        cand_clusters.kind=whole_table_other.kind(cand_clusters_id);
        cand_clusters.product_type=whole_table_other.product(cand_clusters_id); %产物类型可能有多个
        cand_clusters.Accession=whole_table_other.Accession(cand_clusters_id);
        cand_clusters.region_ids=zeros(cand_clusters.num,1); %
    end
    for f=1:length(filesinfolder)
        gbkfilename=filesinfolder(f).name;
        if printprogress
            fprintf('%s\n',gbkfilename)
        end
        clustergbk_path=[loc_folderpath,sepstr,gbkfilename];
        [feature_cds,feature_other,~,~,loc_gbk] = Data_GBK2CDSMisc(clustergbk_path);
        table_cds=struct2table(feature_cds,'AsArray',true);
        table_other=struct2table(feature_other,'AsArray',true);
        CDSs_n=size(table_cds,1);
        for CDS_locs_id = 1:CDSs_n
            if isempty(table_cds.gene_kind{CDS_locs_id})
                table_cds.gene_kind{CDS_locs_id}='other';
            end
        end
        if regions_flag
            % basic information
            regions.num=regions.num+1;
            regions.foldername{regions.num,1}=loc_folderpath;
            regions.ntseq{regions.num,1}=loc_gbk.Sequence;
            regions.seqlen(regions.num,1)=str2double(loc_gbk.LocusSequenceLength);
            regions.gbk_header(regions.num,1)=rmfield(loc_gbk,{'CDS','Sequence','Features'}); % 用更全的gbk_header替代之前的LOCUS_DEF_ACC_VER_strs
            regions.filename{regions.num,1}=gbkfilename(1:(end-4));
            regions.inforstr{regions.num,1}=[loc_folderpath,'_',gbkfilename(1:(end-4)),'_',num2str(loc_gbk.LocusSequenceLength)];
            regions.LocusName{regions.num,1}=loc_gbk.LocusName;
            loc_protocluster_inex=ismember(table_other.Category,'protocluster');
            regions.product_single_type{regions.num,1}=sort(table_other.product(loc_protocluster_inex));
            if contains(lower(loc_gbk.Definition),'plasmid')
                regions.plasmid_flag(regions.num,1)=1;
            else
                regions.plasmid_flag(regions.num,1)=0;
            end
            regions.fulllen(regions.num,1)=wholegbk_LocusSequenceLength(ismember(wholegbk_LocusName,loc_gbk.LocusName));% 需要读整个的gbk，然后再根据regions.LocusName赋值。这个主要是为了后面分析次级代谢在基因组上的位置分布用的。
            regions.LocusTopology(regions.num,1)=wholegbk_LocusTopology(ismember(wholegbk_LocusName,loc_gbk.LocusName)); % 需要读整个的gbk的LocusSequenceLength，然后再根据regions.LocusName赋值。看这个region是在线性还是环形DNA上的
            for i = 1:size(loc_gbk.Comment,1) % if Comment starts with a space line, it will make genbankread can't return comment correctly, you need remove the space line of comment in the gbk file
                if startsWith(loc_gbk.Comment(i,:),'Orig. start  :: ')
                    match = regexp(loc_gbk.Comment(i,:), '\d+', 'match');
                    regions.Orig_start(regions.num,1)=str2double(match{1});
                elseif startsWith(loc_gbk.Comment(i,:),'Orig. end    :: ')
                    match = regexp(loc_gbk.Comment(i,:), '\d+', 'match');
                    regions.Orig_end(regions.num,1)=str2double(match{1});
                end
            end
            loc_region_inex=ismember(table_other.Category,'region');
            if strcmp(table_other.contig_edge{loc_region_inex},'True')
                regions.contig_edge(regions.num,1)=1;
            else
                regions.contig_edge(regions.num,1)=0;
            end
            if cand_clusters_flag
                cand_clusters.region_ids(ismember(cand_clusters.Accession,loc_gbk.Accession)&cand_clusters.border(:,1)>=regions.Orig_start(regions.num,1)&cand_clusters.border(:,2)<=regions.Orig_end(regions.num,1))=regions.num;
            end
            % resistance
            resistance_index=find(strcmp(table_cds.gene_kind, 'resistance'));
            resist=[];
            resist.num=length(resistance_index);
            if resist.num>0
                for resistance_id = 1:length(resistance_index)
                    if ischar(table_cds.gene_functions{resistance_index(resistance_id)})% It maybe contains other gene function annotation
                        resistance_content=table_cds.gene_functions{resistance_index(resistance_id)};
                    else
                        resistance_content=table_cds.gene_functions{resistance_index(resistance_id)}{contains(table_cds.gene_functions{resistance_index(resistance_id)},'resistance (resist)')};
                    end
                    matches = regexp(resistance_content, resist_pattern, 'tokens', 'once');
                    resist.className{resistance_id,1} = matches{1};
                    resist.score(resistance_id,1) = str2double(matches{2});
                    resist.E_value(resistance_id,1) = str2double(matches{3});
                end
            else
                resist.className=[];resist.score=[];resist.E_value=[];
            end
            regions.resist(regions.num,1)=resist;
            % tfbs and tta
            TFBS_struct=[];
            TFBS_struct.num=0;
            TFBS_struct.name=[];TFBS_struct.description=[];TFBS_struct.confidence=[];TFBS_struct.score=[];TFBS_struct.iscomplement=[];TFBS_struct.border=[];TFBS_struct.length=[];TFBS_struct.sequence=[];
            TTA_struct=[];
            TTA_struct.num=0;
            TTA_struct.iscomplement=[];TTA_struct.border=[];TTA_struct.CDSs_id=[];
            misc_feature_id=find(ismember(table_other.Category,'misc_feature'));
            if ~isempty(misc_feature_id)
                misc_feature_note=table_other.note(misc_feature_id);
                misc_feature_note=misc_feature_note(~cellfun('isempty', misc_feature_note));% remove misc_feature which doesn't contain note field
                TFBS_index=find(startsWith(misc_feature_note,'TFBS match to'));
                TFBS_struct.num=length(TFBS_index);
                if TFBS_struct.num>0
                    for i = 1:TFBS_struct.num
                        matches = regexp(misc_feature_note{TFBS_index(i)}, TFBS_pattern, 'names');
                        TFBS_struct.name{i,1}=matches.TFBS;
                        TFBS_struct.description{i,1}=matches.Description;
                        TFBS_struct.confidence{i,1}=matches.Confidence;
                        TFBS_struct.score(i,1)=str2double(matches.Score);
                        if strcmp(table_other.Direction{misc_feature_id(TFBS_index(i))},'reverse')
                            TFBS_struct.iscomplement(i,1)=1;
                        else
                            TFBS_struct.iscomplement(i,1)=0;
                        end
                        TFBS_struct.border(i,[1,2])=[table_other.Minimum(misc_feature_id(TFBS_index(i))),table_other.Maximum(misc_feature_id(TFBS_index(i)))];
                        TFBS_struct.length(i,1)=table_other.Length(misc_feature_id(TFBS_index(i)));
                        if TFBS_struct.iscomplement(i,1)~=1%正向序列，如果是反向的，就互补一下
                            TFBS_struct.sequence{i,1}=loc_gbk.Sequence(TFBS_struct.border(i,1):TFBS_struct.border(i,2));
                        else
                            TFBS_struct.sequence{i,1}=seqrcomplement(loc_gbk.Sequence(TFBS_struct.border(i,1):TFBS_struct.border(i,2)));
                        end
                    end
                end
                TTA_index=find(startsWith(misc_feature_note,'tta leucine codon'));
                TTA_struct.num=length(TTA_index);
                if TTA_struct.num > 0
                    for i = 1:TTA_struct.num
                        if strcmp(table_other.Direction{misc_feature_id(TTA_index(i))},'reverse')
                            TTA_struct.iscomplement(i,1)=1;
                        else
                            TTA_struct.iscomplement(i,1)=0;
                        end
                        TTA_struct.border(i,[1,2])=[table_other.Minimum(misc_feature_id(TTA_index(i))),table_other.Maximum(misc_feature_id(TTA_index(i)))];    
                        for CDS_locs_id = 1:CDSs_n
                            if table_cds.Minimum(CDS_locs_id)<=TTA_struct.border(i,1)&&TTA_struct.border(i,2)<=table_cds.Maximum(CDS_locs_id)
                               TTA_struct.CDSs_id(i,1)=CDS_locs_id;
                               break
                            end
                        end
                    end
                end
            end
            regions.TFBS(regions.num,1)=TFBS_struct;
            regions.TTA(regions.num,1)=TTA_struct;
            % MiBiG comparison
            jsonData_id=0;
            loc_region_locus_name=split(regions.LocusName{regions.num},'.');
            loc_region_locus_name=loc_region_locus_name{1};
            for i = 1:length(jsonData.records)
                loc_json_locus_name=split(jsonData.records(i).name,'.');
                loc_json_locus_name=loc_json_locus_name{1};
                if strcmp(loc_json_locus_name,loc_region_locus_name)
                    jsonData_id=i;
                    break
                end
            end
            matches = regexp(regions.filename{regions.num,1}, region_pattern, 'tokens');
            region_id_in_json=matches{1}{1};%str
            if isfield(jsonData.records(jsonData_id).modules,'antismash_modules_cluster_compare')
                % ProtoToRegion_RiQ=jsonData.records(jsonData_id).modules.antismash_modules_cluster_compare.db_results.MIBiG.by_region.(['x',region_id_in_json]).ProtoToRegion_RiQ;
                ProtoToRegion_RiQ=jsonData.records(jsonData_id).modules.antismash_modules_cluster_compare.db_results.(comparison_data).by_region.(['x',region_id_in_json]).ProtoToRegion_RiQ; % use custom_db for comparison
                mibig_comparison=[];
                mibig_id_raw=fieldnames(ProtoToRegion_RiQ.scores_by_region);
                mibig_comparison.num=length(mibig_id_raw);
                if mibig_comparison.num>0
                    for i = 1:mibig_comparison.num
                        BGC_raw=split(mibig_id_raw{i},'_');
                        mibig_comparison.mibig_id{i,1}=BGC_raw{1};
                        mibig_comparison.border(i,:)=str2double(BGC_raw(2:3));
                        mibig_comparison.similarity_score(i,1)=ProtoToRegion_RiQ.scores_by_region.(mibig_id_raw{i});
                        if strcmp(comparison_data,'MIBiG')
                            mibig_comparison.compound{i,1}=ProtoToRegion_RiQ.reference_regions.(mibig_id_raw{i}).description;
                            mibig_comparison.product_type{i,1}=ProtoToRegion_RiQ.reference_regions.(mibig_id_raw{i}).products;
                            mibig_comparison.organism{i,1}=ProtoToRegion_RiQ.reference_regions.(mibig_id_raw{i}).organism;
                        end
                    end
                else
                    mibig_comparison.mibig_id=[];mibig_comparison.border=[];mibig_comparison.similarity_score=[];
                    if strcmp(comparison_data,'MIBiG')
                        mibig_comparison.compound=[];mibig_comparison.product_type=[];mibig_comparison.organism=[];
                    end
                end
                regions.mibig_comparison(regions.num,1)=mibig_comparison;
            end
            % knownclusterblast
            if isfield(jsonData.records(jsonData_id).modules,'antismash_modules_clusterblast')
                region_id_in_knownclusterblast=0;
                blast_results_all=jsonData.records(jsonData_id).modules.antismash_modules_clusterblast.knowncluster.results;
                for i = 1:length(blast_results_all)
                    if blast_results_all(i).region_number==str2double(region_id_in_json)
                        region_id_in_knownclusterblast=i;
                        break
                    end
                end
                blast_results=blast_results_all(region_id_in_knownclusterblast);
                knownclusterblast=[];
                knownclusterblast.num=blast_results.total_hits;
                if knownclusterblast.num>0
                    for i = 1:length(blast_results.ranking) %blast_results.ranking最多只有50个，即使knownclusterblast.num>50,所以用length(blast_results.ranking)而不是knownclusterblast.num
                        knownclusterblast.mibig_id{i,1}=blast_results.ranking{i}{1}.accession;
                        knownclusterblast.compound{i,1}=blast_results.ranking{i}{1}.description;
                        knownclusterblast.product_type{i,1}=blast_results.ranking{i}{1}.cluster_type; %多种type是用+连接的
                        knownclusterblast.blast_score(i,1)=blast_results.ranking{i}{2}.blast_score;
                    end
                else
                    knownclusterblast.mibig_id=[];knownclusterblast.compound=[];knownclusterblast.product_type=[];knownclusterblast.blast_score=[];
                end
                regions.knownclusterblast(regions.num,1)=knownclusterblast;
            end
            %t2pks
            protocluster_id=find(ismember(table_other.Category,'protocluster'));
            t2pks_id=protocluster_id(ismember(table_other.product(protocluster_id),'T2PKS'));
            t2pks=[];
            t2pks.num=length(t2pks_id);
            if t2pks.num>0
                for i = 1:t2pks.num
                    t2pks.border(i,:)=[table_other.Minimum(t2pks_id(i)),table_other.Maximum(t2pks_id(i))];
                    if strcmp(table_other.Direction{t2pks_id(i)},'reverse')
                        t2pks.iscomplement(i,1)=1;
                    else
                        t2pks.iscomplement(i,1)=0;
                    end
                    if strcmp(table_other.contig_edge{t2pks_id(i)},'True')
                        t2pks.contig_edge(i,1)=1;
                    else
                        t2pks.contig_edge(i,1)=0;
                    end
                    if ismember('t2pks_malonyl_elongations',table_other.Properties.VariableNames) % not everyone hase this information
                        if ischar(table_other.t2pks_malonyl_elongations{t2pks_id(i)})
                            molecule_id=1;
                            match = regexp(table_other.t2pks_malonyl_elongations{t2pks_id(i)}, t2pks_elongation_pattern, 'tokens', 'once');
                            t2pks.malonyl_elongations_n{i,1}{molecule_id,1}=match{1};
                            t2pks.malonyl_elongations_score{i,1}{molecule_id,1}=str2double(match{2});% need change from double to cell
                            t2pks.malonyl_elongations_E_value{i,1}{molecule_id,1}=str2double(match{3});
                        else
                            molecule_n=length(table_other.t2pks_malonyl_elongations{t2pks_id(i)});
                            for molecule_id = 1:molecule_n
                                match = regexp(table_other.t2pks_malonyl_elongations{t2pks_id(i)}{molecule_id}, t2pks_elongation_pattern, 'tokens', 'once');
                                t2pks.malonyl_elongations_n{i,1}{molecule_id,1}=match{1};
                                t2pks.malonyl_elongations_score{i,1}{molecule_id,1}=str2double(match{2});
                                t2pks.malonyl_elongations_E_value{i,1}{molecule_id,1}=str2double(match{3});
                            end
                        end
                    else
                        t2pks.malonyl_elongations_n{i,1}=[];
                        t2pks.malonyl_elongations_score{i,1}=[];
                        t2pks.malonyl_elongations_E_value{i,1}=[];
                    end
                    if ismember('t2pks_molecular_weights',table_other.Properties.VariableNames) % not everyone hase this information
                        if ischar(table_other.t2pks_molecular_weights{t2pks_id(i)}) % if only one result, its type will be char, not cell
                            molecule_id=1;
                            loc_molecule=split(table_other.t2pks_molecular_weights{t2pks_id(i)},': ');
                            t2pks.molecular_weights{i,1}{molecule_id,1}=loc_molecule{1};
                            t2pks.molecular_weights{i,1}{molecule_id,2}=str2double(loc_molecule{2});
                        else
                            molecule_n=length(table_other.t2pks_molecular_weights{t2pks_id(i)});
                            for molecule_id = 1:molecule_n
                                loc_molecule=split(table_other.t2pks_molecular_weights{t2pks_id(i)}{1},': ');
                                t2pks.molecular_weights{i,1}{molecule_id,1}=loc_molecule{1};
                                t2pks.molecular_weights{i,1}{molecule_id,2}=str2double(loc_molecule{2});
                            end
                        end
                    else
                        t2pks.molecular_weights{i,1}=[];
                    end
                    if ismember('t2pks_product_classes',table_other.Properties.VariableNames) % not everyone hase this information
                        t2pks.product_classes{i,1}=table_other.t2pks_product_classes{t2pks_id(i)};
                    else
                        t2pks.product_classes{i,1}=[];
                    end
                    if ismember('t2pks_starter_units',table_other.Properties.VariableNames) % not everyone hase this information
                        if ischar(table_other.t2pks_starter_units{t2pks_id(i)}) % if only one result, its type will be char, not cell
                            molecule_id=1; % it could be more than one
                            match = regexp(table_other.t2pks_starter_units{t2pks_id(i)}, t2pks_starter_units_pattern, 'tokens', 'once');
                            t2pks.starter_units{i,1}{molecule_id,1}=match{1};
                            t2pks.starter_units_score{i,1}(molecule_id,1)=str2double(match{2});
                            t2pks.starter_units_E_value{i,1}(molecule_id,1)=str2double(match{3});
                        else
                            molecule_n=length(table_other.t2pks_starter_units{t2pks_id(i)}); % it could be more than one
                            for molecule_id = 1:molecule_n
                                match = regexp(table_other.t2pks_starter_units{t2pks_id(i)}{molecule_id}, t2pks_starter_units_pattern, 'tokens', 'once');
                                t2pks.starter_units{i,1}{molecule_id,1}=match{1};
                                t2pks.starter_units_score{i,1}(molecule_id,1)=str2double(match{2});
                                t2pks.starter_units_E_value{i,1}(molecule_id,1)=str2double(match{3});
                            end
                        end
                    else
                        t2pks.starter_units{i,1}=[];
                        t2pks.starter_units_score{i,1}=[];
                        t2pks.starter_units_E_value{i,1}=[];
                    end
                end
            else
                t2pks.border=[];t2pks.iscomplement=[];t2pks.contig_edge=[];t2pks.malonyl_elongations_n=[];t2pks.malonyl_elongations_score=[];t2pks.malonyl_elongations_E_value=[];t2pks.molecular_weights=[];t2pks.product_classes=[];t2pks.starter_units=[];t2pks.starter_units_score=[];t2pks.starter_units_E_value=[];
            end
            regions.t2pks(regions.num,1)=t2pks;
            % RiPP 只是antismash预测有核心肽的RiPP，并不是所有的RiPP都有记录
            ripp=[]; %记录lanthipeptides,lassopeptides,sactipeptides,thiopeptides四类RiPP的肽段。有的虽然是这四类RiPP，但是没有预测的肽段。注意RiPP不止这4种。
            ripp.num=0;
            ripp.infor=[];
            if ~isempty(ripp_all_list{1})
                if ismember('locus_tag',table_cds.Properties.VariableNames)
                    loc_locus_tag=table_cds.locus_tag;
                    loc_locus_tag(cellfun(@isempty,loc_locus_tag))=[];
                    loc_ripp_index_by_locus_tag=ismember(ripp_all_list{1},loc_locus_tag);
                else
                    loc_ripp_index_by_locus_tag=zeros(length(ripp_all_list{1}),1);
                end
                if ismember('gene',table_cds.Properties.VariableNames)
                    loc_gene=table_cds.gene;
                    loc_gene(cellfun(@isempty,loc_gene))=[];
                    loc_ripp_index_by_gene=ismember(ripp_all_list{1},loc_gene);
                else
                    loc_ripp_index_by_gene=zeros(length(ripp_all_list{1}),1);
                end
                if ismember('protein_id',table_cds.Properties.VariableNames)
                    loc_protein_id=table_cds.protein_id;
                    loc_protein_id(cellfun(@isempty,loc_protein_id))=[];
                    loc_ripp_index_by_protein_id=ismember(ripp_all_list{1},loc_protein_id);
                else
                    loc_ripp_index_by_protein_id=zeros(length(ripp_all_list{1}),1);
                end
                if ismember('locus_tag',table_cds.Properties.VariableNames)||ismember('gene',table_cds.Properties.VariableNames)||ismember('protein_id',table_cds.Properties.VariableNames)
                    loc_ripp_index_by_gene=reshape(loc_ripp_index_by_gene,size(loc_ripp_index_by_locus_tag));
                    loc_ripp_index_by_protein_id=reshape(loc_ripp_index_by_protein_id,size(loc_ripp_index_by_locus_tag));
                    loc_ripp_index=loc_ripp_index_by_locus_tag|loc_ripp_index_by_gene|loc_ripp_index_by_protein_id;
                else
                    error('Check which field is used in RiPP information!')
                end
                if any(loc_ripp_index)
                    ripp.num=sum(loc_ripp_index);
                    ripp.infor=ripp_all_list{2}(loc_ripp_index);
                end
            end
            regions.ripp_with_infor(regions.num,1)=ripp;
            % CompaRiPPson的信息没有提
            % merge candidate cluster
            cand_cluster_tmp=[];
            loc_index=find(ismember(table_other.Category,'cand_cluster'));
            loc_protocluster_num=sum(ismember(table_other.Category,'protocluster'));
            cand_cluster_tmp.kind=table_other.kind(loc_index);% 第一次，给判定cand_cluster的kind
            hybrid_num=sum(ismember(cand_cluster_tmp.kind,'chemical_hybrid'));% 2 protocluster in chemical_hybrid
            interleaved_num=sum(ismember(cand_cluster_tmp.kind,'interleaved'));% 2 protocluster in interleaved
            if length(loc_index)>1&&length(loc_index)~=(loc_protocluster_num-hybrid_num-interleaved_num)%In general, there is no neighbouring type. But if the only one type is neighbouring, save it.
                loc_index=loc_index(~ismember(cand_cluster_tmp.kind,'neighbouring'));
                cand_cluster_tmp.kind=table_other.kind(loc_index);% 如果loc_index可能有改动，需要第二次给cand_cluster的kind赋值
            end
            % cand_cluster_tmp.cand_cluster_ids=str2double(table_other.candidate_cluster_number(loc_index));
            cand_cluster_tmp.product=table_other.product(loc_index);
            cand_cluster_tmp.protoclusters=table_other.protoclusters(loc_index);% in the region
            cand_cluster_tmp.Minimum=table_other.Minimum(loc_index);
            cand_cluster_tmp.Maximum=table_other.Maximum(loc_index);
            initial_index=1;
            for i = 1:length(loc_index)
                protoclusters_index=find(cand_clusters.region_ids==regions.num&ismember(cand_clusters.border,[cand_cluster_tmp.Minimum(i)+regions.Orig_start(regions.num,1),cand_cluster_tmp.Maximum(i)+regions.Orig_start(regions.num,1)],'rows'));
                cand_cluster_tmp.cand_cluster_ids(i)=cand_clusters.candidate_cluster_number(protoclusters_index(initial_index));
                cand_cluster_tmp.protoclusters_index(i)=protoclusters_index(initial_index);% only use in merge product for finding protocluster id in mibig comparison
                if length(protoclusters_index)>1
                    initial_index=initial_index+1; % It's possible that there are more than one protoclusters in the same location and region
                end
                if initial_index>length(protoclusters_index) % If excess the range, restart from 1
                    initial_index=1;
                end
            end
            if length(loc_index)>1 % only analyze region with more than one products
                protocluster=[];
                protocluster=extract_infor(protocluster,table_other,'protocluster');
                proto_core=[];
                proto_core=extract_infor(proto_core,table_other,'proto_core');
                loc_pair=[];
                loc_cand_id_pair=[];
                loc_neighbourhood_overlap_ratio=[];
                for i = 1:length(cand_cluster_tmp.cand_cluster_ids)
                    for k = i+1:length(cand_cluster_tmp.cand_cluster_ids)
                        loc_pair=[loc_pair;[i,k]];
                        loc_cand_id_pair=[loc_cand_id_pair;[cand_cluster_tmp.cand_cluster_ids(i),cand_cluster_tmp.cand_cluster_ids(k)]];
                        loc_neighbourhood_overlap_ratio=[loc_neighbourhood_overlap_ratio;max([cal_max_overlap_ratio(cand_cluster_tmp,proto_core,k,i,protocluster),cal_max_overlap_ratio(cand_cluster_tmp,proto_core,i,k,protocluster)])];% ratio的时候选大的那个。A对B，不等于B对A。
                    end
                end
                close_pair=loc_pair(loc_neighbourhood_overlap_ratio>=merge_threshold,:);
                far_pair=loc_pair(loc_neighbourhood_overlap_ratio<merge_threshold,:);
                merge_result_id=1;
                merge_result=[];
                while ~isempty(close_pair)%把近的进行合并
                    merge_result{merge_result_id}=close_pair(1,:);
                    close_pair(1,:)=[];
                    remove_list=[];
                    for i=1:size(close_pair,1)
                        if ~isempty(intersect(merge_result{merge_result_id},close_pair(i,:)))
                            merge_result{merge_result_id}=union(merge_result{merge_result_id},close_pair(i,:));
                            remove_list=[remove_list;i];
                        end
                    end
                    close_pair(remove_list,:)=[];
                    merge_result_id=merge_result_id+1;
                end
                single_result=[];
                for i = 1:size(far_pair,1)%记录距离远的id
                    single_result=union(single_result,far_pair(i,:));
                end
                for i = 1:length(merge_result)%去掉在合并的里面已经出现过的
                    single_result=setdiff(single_result,merge_result{i});
                end
            else
                merge_result=[];
                single_result=1;
            end
            product_after_merge=[];
            product_after_merge.total_num=length(merge_result)+length(single_result);
            product_after_merge.merge_num=length(merge_result);
            product_after_merge.unmerge_num=length(single_result);
            product_after_merge.merge_ratio=1-product_after_merge.unmerge_num/length(cand_cluster_tmp.product);% merged candidate clusters (without neighbouring)/totoal candidate clusters (without neighbouring). In general, there is no neighbouring type. But if the only one type is neighbouring, save it.
            product_after_merge.ismerge=[ones(product_after_merge.merge_num,1);zeros(product_after_merge.unmerge_num,1)];
            for i = 1:product_after_merge.merge_num%现在这里是在region上的位置，在整个基因组上的位置需要加上region.Orig_start
                product_after_merge.border(i,:)=[min(cand_cluster_tmp.Minimum(merge_result{i})),max(cand_cluster_tmp.Maximum(merge_result{i}))];
                product_after_merge.cand_cluster_ids{i,1}=cand_cluster_tmp.cand_cluster_ids(merge_result{i});% cand_clusters中的id
                product_after_merge.cand_cluster_kind{i,1}=cand_cluster_tmp.kind(merge_result{i});
                product_after_merge.protocluster_ids{i,1}=[];% protocluster_ids是当前region内部的protocluster序号
                product_after_merge.protoclusters_index{i,1}=[];% 当一个文件夹里有多个contig时，protoclusters_index用于找到整个大gbk(foldername.gbk)中正确的global protoclusters id, only use in merge product for finding protocluster id in mibig comparison
                product_after_merge.cand_cluster_product_type{i,1}=[]; % (与product_type相比) cand_cluster_product是merge之前的，一个cand_cluster放一个cell，一个cand_cluster可能有多种product_type
                product_after_merge.product_type{i,1}=[]; %(与cand_cluster_product相比) product_type是merge以后的，相当于把cand_cluster_product里面的每个product都放在一个cell里了
                for k = 1:length(merge_result{i})
                    product_after_merge.protocluster_ids{i}=[product_after_merge.protocluster_ids{i};str2double(cand_cluster_tmp.protoclusters{merge_result{i}(k)})];
                    product_after_merge.protoclusters_index{i}=[product_after_merge.protoclusters_index{i};cand_cluster_tmp.protoclusters_index(merge_result{i}(k))];
                    if ischar(cand_cluster_tmp.product{merge_result{i}(k)})
                        product_after_merge.cand_cluster_product_type{i,1}{k,1}=cand_cluster_tmp.product(merge_result{i}(k));
                        product_after_merge.product_type{i,1}=[product_after_merge.product_type{i};cand_cluster_tmp.product(merge_result{i}(k))];
                    else
                        product_after_merge.cand_cluster_product_type{i,1}{k,1}=cand_cluster_tmp.product{merge_result{i}(k)};
                        product_after_merge.product_type{i,1}=[product_after_merge.product_type{i};cand_cluster_tmp.product{merge_result{i}(k)}];
                    end
                end
                product_after_merge.protocluster_ids{i}=unique(product_after_merge.protocluster_ids{i}); % due to some neighbouring kind candidate cluster, there may be some repeats
            end
            for i = 1:product_after_merge.unmerge_num
                single_index=i+product_after_merge.merge_num;
                product_after_merge.border(single_index,:)=[min(cand_cluster_tmp.Minimum(single_result(i))),max(cand_cluster_tmp.Maximum(single_result(i)))];
                product_after_merge.cand_cluster_ids{single_index,1}=cand_cluster_tmp.cand_cluster_ids(single_result(i));% cand_clusters中的id
                product_after_merge.cand_cluster_kind{single_index,1}=cand_cluster_tmp.kind(single_result(i));
                product_after_merge.protocluster_ids{single_index,1}=str2double(cand_cluster_tmp.protoclusters{single_result(i)});
                product_after_merge.protoclusters_index{single_index,1}=cand_cluster_tmp.protoclusters_index(single_result(i));
                if ischar(cand_cluster_tmp.product{single_result(i)})
                    product_after_merge.cand_cluster_product_type{single_index,1}=cand_cluster_tmp.product(single_result(i));
                    product_after_merge.product_type{single_index,1}=cand_cluster_tmp.product(single_result(i));% 虽然没有经过合并，但仍然可能有多个产物，如chemical hybird
                else
                    product_after_merge.cand_cluster_product_type{single_index,1}=cand_cluster_tmp.product{single_result(i)};
                    product_after_merge.product_type{single_index,1}=cand_cluster_tmp.product{single_result(i)};% 虽然没有经过合并，但仍然可能有多个产物，如chemical hybird
                end
            end
            for i = 1:product_after_merge.total_num
                uni_sorted_product=unique(product_after_merge.product_type{i},'sorted');
                product_after_merge.uni_sorted_product_type_str(i,1)=join(uni_sorted_product,'+');
            end
            % reorder
            if product_after_merge.total_num>1
                [~,I]=sort(product_after_merge.border(:,1),'ascend');
                product_after_merge.border=product_after_merge.border(I,:);
                product_after_merge.cand_cluster_ids=product_after_merge.cand_cluster_ids(I);
                product_after_merge.cand_cluster_kind=product_after_merge.cand_cluster_kind(I);
                product_after_merge.protocluster_ids=product_after_merge.protocluster_ids(I);
                product_after_merge.protoclusters_index=product_after_merge.protoclusters_index(I);
                product_after_merge.cand_cluster_product_type=product_after_merge.cand_cluster_product_type(I);
                product_after_merge.product_type=product_after_merge.product_type(I);
                product_after_merge.uni_sorted_product_type_str=product_after_merge.uni_sorted_product_type_str(I);
                % unmerged_max_neighbourhood_overlap_ratio用于check不合并的组合之间最大的邻域重叠率的分布
                for i = 1:product_after_merge.total_num-1
                    cand_id1=product_after_merge.cand_cluster_ids{i};
                    cand_id2=product_after_merge.cand_cluster_ids{i+1};
                    loc_neighbourhood_overlap_ratio_list=[];
                    for cand_id1_index = 1:length(cand_id1)
                        for cand_id2_index = 1:length(cand_id2)
                            loc_neighbourhood_overlap_ratio_list=[loc_neighbourhood_overlap_ratio_list;loc_neighbourhood_overlap_ratio(ismember(loc_cand_id_pair,[cand_id1(cand_id1_index),cand_id2(cand_id2_index)],'rows'))];
                        end
                    end
                    product_after_merge.unmerged_max_neighbourhood_overlap_ratio(i,1)=max(loc_neighbourhood_overlap_ratio_list);% length=product_after_merge.total_num-1
                    middle_border=mean([product_after_merge.border(i,2),product_after_merge.border(i+1,1)]);
                    for CDS_locs_id = 1:CDSs_n
                        if table_cds.Minimum(CDS_locs_id)<=middle_border&&middle_border<=table_cds.Maximum(CDS_locs_id)
                            if (middle_border-table_cds.Minimum(CDS_locs_id))/table_cds.Length(CDS_locs_id)>=0.5
                                product_after_merge.border(i,2)=table_cds.Maximum(CDS_locs_id);
                                product_after_merge.border(i+1,1)=table_cds.Maximum(CDS_locs_id)+1;
                            else
                                product_after_merge.border(i,2)=table_cds.Minimum(CDS_locs_id)-1;
                                product_after_merge.border(i+1,1)=table_cds.Minimum(CDS_locs_id);
                            end
                        break
                        end
                    end
                end
            else
                product_after_merge.unmerged_max_neighbourhood_overlap_ratio=[];
            end
            product_after_merge.border(1,1)=table_cds.Minimum(find(product_after_merge.border(1,1)<=table_cds.Minimum,1,'first'));
            product_after_merge.border(end,2)=table_cds.Maximum(find(product_after_merge.border(end,2)>=table_cds.Maximum,1,'last'));
            for i = 1:product_after_merge.total_num
                % product
                if isfield(jsonData.records(jsonData_id).modules,'antismash_modules_cluster_compare')
                    product_after_merge.protoclusters_in_single_seq{i,1}=[];
                    for k = 1:length(product_after_merge.cand_cluster_ids{i})
                         product_after_merge.protoclusters_in_single_seq{i,1}=[product_after_merge.protoclusters_in_single_seq{i,1};cand_clusters.protoclusters{product_after_merge.protoclusters_index{i}(k)}];% if use protoclusters, the order is wrong when there more than one contig sequence
                    end
                    product_after_merge.protoclusters_in_single_seq{i,1}=unique(product_after_merge.protoclusters_in_single_seq{i,1}); % due to some neighbouring kind candidate cluster, there may be some repeats
                    product_after_merge.protocluster_without_mibig_ratio(i,1)=0;% 该merged product含有没有mibig结果的protocluster
                    product_after_merge.bad_merge_ratio(i,1)=0;
                    blank_flag=0;
                    if ~isempty(fieldnames(ProtoToRegion_RiQ.scores_by_region))
                        proto_num=length(product_after_merge.protoclusters_in_single_seq{i,1});
                        mibig_infor_list=cell(proto_num,1);%1.id 2.border 3.identity 4.order 5.component 6.final_score
                        for k = 1:proto_num
                            loc_field=['x',num2str(product_after_merge.protoclusters_in_single_seq{i,1}(k))];
                            if isfield(ProtoToRegion_RiQ.details.details,loc_field)% 有可能region里有mibig comparison结果，但是某个protoclust没有mibig comparison结果
                                proto_mibig=ProtoToRegion_RiQ.details.details.(loc_field);
                                proto_mibig_id_raw=fieldnames(proto_mibig);
                                proto_mibig_comparison_num=length(proto_mibig_id_raw);
                                if proto_mibig_comparison_num>0
                                    for loc_proto_mibig_index = 1:proto_mibig_comparison_num
                                        proto_BGC_raw=split(proto_mibig_id_raw{loc_proto_mibig_index},'_');
                                        mibig_infor_list{k}{1}(loc_proto_mibig_index,1)=join(proto_BGC_raw(1:end-2),'_');
                                        mibig_infor_list{k}{2}(loc_proto_mibig_index,1)=proto_mibig.(proto_mibig_id_raw{loc_proto_mibig_index}).final_score;
                                    end
                                else
                                    mibig_infor_list{k}{1}=[];mibig_infor_list{k}{2}=[];
                                end
                            else
                                mibig_infor_list{k}{1}=[];mibig_infor_list{k}{2}=[];
                            end
                        end
                        if proto_num~=1||~isempty(mibig_infor_list{1}{1})%如果只有1个protocluster，且没有mibig comparison结果，就留空
                            mibig_id_final=[];
                            % 取并集，记录没有mibig comparison结果的protocluster，并进行标注product_after_merge.protocluster_without_mibig_ratio(i,1);其值在0~1之间
                            % 差异越大的protocluster，最后并集的元素数量比单个的最大值的数量会越多。这一项用多出来的元素个数除以单个protocluster最多的元素个数表示
                            % 并进行标注product_after_merge.bad_merge_ratio(i,1)。值为非负数，可能会大于1
                            protocluster_without_mibig_n=0;
                            for k = 1:proto_num
                                if ~isempty(mibig_infor_list{k}{1}) % 只使用具有mibig comparison的protocluster的final score等信息进行计算
                                    mibig_id_final=union(mibig_id_final,mibig_infor_list{k}{1});
                                else
                                    protocluster_without_mibig_n=protocluster_without_mibig_n+1;
                                end
                            end
                            if proto_num-protocluster_without_mibig_n>0 % 有含有mibig comprison结果的protocluster
                                mibig_infor_num=zeros(proto_num,1);
                                mibig_infor=cell(proto_num,1);
                                for k = 1:proto_num
                                    mibig_infor_num(k)=length(mibig_infor_list{k}{1});
                                    mibig_infor{k}=mibig_infor_list{k}{1};
                                end
                                product_after_merge.protocluster_without_mibig_ratio(i,1)=protocluster_without_mibig_n/proto_num;
                                product_after_merge.bad_merge_ratio(i,1)=length(mibig_id_final)/max(mibig_infor_num);
                                mibig_score_final_modified=zeros(length(mibig_id_final),2);
                                mibig_score_final_raw=zeros(length(mibig_id_final),proto_num);
                                for k = 1:length(mibig_id_final)
                                    for loc_proto_mibig_index = 1:proto_num
                                        loc_mibig_index=ismember(mibig_infor_list{loc_proto_mibig_index}{1},mibig_id_final{k});
                                        if any(loc_mibig_index)
                                            mibig_score_final_raw(k,loc_proto_mibig_index)=mibig_infor_list{loc_proto_mibig_index}{2}(loc_mibig_index);
                                        end
                                        mibig_score_final_modified(k,2)=sum(exp(3*mibig_score_final_raw(k,:))-1)/(exp(3)-1)/proto_num; % similar score normalized by protocluster number
                                    end
                                end
                                mibig_score_final_modified(:,1)=mean(mibig_score_final_raw,2); % simple average
                                [~,I]=sort(mibig_score_final_modified(:,2),'descend');
                                mibig_id_final=mibig_id_final(I);
                                product_after_merge.mibig_id_final{i,1}=mibig_id_final;
                                product_after_merge.mibig_score_final_raw{i,1}=mibig_score_final_raw(I,:);
                                product_after_merge.mibig_score_final_modified{i,1}=mibig_score_final_modified(I,:);
                                % for custom database no these information
                                if strcmp(comparison_data,'MIBiG')
                                    for k = 1:length(mibig_id_final) % add mibig information 后续用其他数据库比对的话，这里就要改掉了，需要提前输入一个库的对照数据
                                        loc_mibig_index=ismember(regions.mibig_comparison(regions.num).mibig_id,mibig_id_final{k});
                                        if any(loc_mibig_index)
                                            product_after_merge.mibig_compound{i,1}{k,1}=regions.mibig_comparison(regions.num).compound{loc_mibig_index};
                                            product_after_merge.mibig_product_type{i,1}{k,1}=regions.mibig_comparison(regions.num).product_type{loc_mibig_index};
                                            product_after_merge.mibig_organism{i,1}{k,1}=regions.mibig_comparison(regions.num).organism{loc_mibig_index};
                                        else
                                            error('No mibig information')
                                        end
                                    end
                                end
                            else
                                blank_flag=1;
                            end
                        else
                            blank_flag=1;
                        end
                    else
                        blank_flag=1;
                    end
                    if blank_flag==1
                        product_after_merge.protocluster_without_mibig_ratio(i,1)=1;
                        product_after_merge.mibig_id_final{i,1}=[];
                        product_after_merge.mibig_score_final_raw{i,1}=[];
                        product_after_merge.mibig_score_final_modified{i,1}=[];
                        if strcmp(comparison_data,'MIBiG')
                            product_after_merge.mibig_compound{i,1}=[];
                            product_after_merge.mibig_product_type{i,1}=[];
                            product_after_merge.mibig_organism{i,1}=[];
                        end
                    end
                end
                % calculate length
                CDS_locs_id_list=find(product_after_merge.border(i,1)<=table_cds.Minimum,1,'first'):find(product_after_merge.border(i,2)>=table_cds.Maximum,1,'last');
                [product_after_merge.cluster_length_biosynthetic_core(i,1),product_after_merge.cds_length_biosynthetic_core(i,1)] = calculate_BGC_length(table_cds,CDS_locs_id_list,'biosynthetic');% core biosynthetic
                [product_after_merge.cluster_length_biosynthetic_all(i,1),product_after_merge.cds_length_biosynthetic_all(i,1)] = calculate_BGC_length(table_cds,CDS_locs_id_list,'biosynthetic+');% core biosynthetic and additional biosynthetic
                [product_after_merge.cluster_length_functional(i,1),product_after_merge.cds_length_functional(i,1)] = calculate_BGC_length(table_cds,CDS_locs_id_list,'functional');% without "other" kind gene
                [product_after_merge.cluster_length_all(i,1),product_after_merge.cds_length_all(i,1)] = calculate_BGC_length(table_cds,CDS_locs_id_list,'all');% all including "other" kind gene
            end
            regions.product_after_merge(regions.num,1)=product_after_merge;
            regions.cluster_length_biosynthetic_core_sum(regions.num,1)=sum(product_after_merge.cluster_length_biosynthetic_core);
            regions.cluster_length_biosynthetic_all_sum(regions.num,1)=sum(product_after_merge.cluster_length_biosynthetic_all);
            regions.cluster_length_functional_sum(regions.num,1)=sum(product_after_merge.cluster_length_functional);
            regions.cluster_length_all_sum(regions.num,1)=sum(product_after_merge.cluster_length_all);
            regions.cds_length_biosynthetic_core_sum(regions.num,1)=sum(product_after_merge.cds_length_biosynthetic_core);
            regions.cds_length_biosynthetic_all_sum(regions.num,1)=sum(product_after_merge.cds_length_biosynthetic_all);
            regions.cds_length_functional_sum(regions.num,1)=sum(product_after_merge.cds_length_functional);
            regions.cds_length_all_sum(regions.num,1)=sum(product_after_merge.cds_length_all);
        end
        if omains_flag
            asdomain_locs=find(ismember(table_other.Category,'aSDomain'));
            domains_iscomplement=zeros(length(asdomain_locs),1);
            domains_sematrix=cell(length(asdomain_locs),1);
            domains_border=zeros(length(asdomain_locs),2);
            domains_seposis_sepnum=zeros(length(asdomain_locs),1);
            domains_typestr=cell(length(asdomain_locs),1);
            domains_typeid=zeros(length(asdomain_locs),1);
            domains_subtypestr=cell(length(asdomain_locs),1);
            domains_locustag=cell(length(asdomain_locs),1);
            domains_translation=cell(length(asdomain_locs),1);
            domains_ntseqlen=zeros(length(asdomain_locs),1);
            domains_maxoverlapp_cdsid=zeros(length(asdomain_locs),2);
            domains_specificity_mat=cell(length(asdomain_locs),1); % CDS依赖
            for asdomain_locs_id = 1:length(asdomain_locs)
                [domains_locustag{asdomain_locs_id},domains_translation{asdomain_locs_id},domains_iscomplement(asdomain_locs_id),domains_border(asdomain_locs_id,:),domains_seposis_sepnum(asdomain_locs_id),domains_sematrix{asdomain_locs_id},domains_ntseqlen(asdomain_locs_id)] = infor_generate(table_other,asdomain_locs(asdomain_locs_id));
                if isfield(feature_other,'domain_subtype')
                    domains_subtypestr{asdomain_locs_id}=table_other.domain_subtype{asdomain_locs(asdomain_locs_id)};
                else
                    domains_subtypestr{asdomain_locs_id}='';
                end
                domains_typestr{asdomain_locs_id}=table_other.aSDomain{asdomain_locs(asdomain_locs_id)};
                domains_typeid(asdomain_locs_id)=domainnames.dtype((strcmp(domainnames.formal,domains_typestr{asdomain_locs_id})));
                if isfield(feature_other,'specificity')
                    spestr=table_other.specificity{asdomain_locs(asdomain_locs_id)};
                    specimat=cell(size(spestr,1),2);
                    for s = 1:size(spestr,1)
                        if ischar(spestr)
                            spesplitedbycomma=strsplit(spestr,':');
                        else
                            spesplitedbycomma=strsplit(spestr{s},':');
                        end
                        specimat{s,1}=strip(spesplitedbycomma{1});
                        if length(spesplitedbycomma)>1
                            specimat{s,2}=strip(spesplitedbycomma{2});
                        else
                            specimat{s,2}='';
                        end
                    end
                else
                    specimat={'',''};
                end
                domains_specificity_mat{asdomain_locs_id}=specimat;
            end
        end
        if CDSs_flag
            for CDS_locs_id = 1:CDSs_n
                CDSs.num=CDSs.num+1;
                [CDS_locustr,CDS_translation,CDSiscomplement,CDS_border,seposis_sepnum,seposis_matrix,CDS_length] = infor_generate(table_cds,CDS_locs_id);
                if CDS_translation(end)~='*'
                    CDS_translation(end+1)='*';
                end
                CDS_inforstr=[regions.inforstr{regions.num,1},'_',CDS_locustr,'_',num2str(CDSs.num),'_',num2str(length(CDS_translation))];
                aaseqlen=length(CDS_translation);
                % determin the way of translation in this CDS, three varibles are all in the same order, from translation start to end
                id_in_regiontseq=[]; % id in the sequence
                seq_in_transorder='';% the order of translation
                id_in_CDStrans=reshape(repmat(1:aaseqlen,3,1),1,aaseqlen*3); % this two must be the same size
                ntseq_gbkfile=regions.ntseq{regions.num,1};
                seposis_list_mat=seposis_matrix;
                for k=1:size(seposis_list_mat,1)
                    selist_on_nt=seposis_list_mat(k,1):seposis_list_mat(k,2);
                    if ~isempty(ntseq_gbkfile)
                        if CDSiscomplement
                            seq_in_transorder=[seqrcomplement(ntseq_gbkfile(selist_on_nt)),seq_in_transorder];
                            selist_on_nt=seposis_list_mat(k,2):-1:seposis_list_mat(k,1);
                            id_in_regiontseq=[selist_on_nt,id_in_regiontseq];
                        else
                            seq_in_transorder=[seq_in_transorder,ntseq_gbkfile(selist_on_nt)];
                            id_in_regiontseq=[id_in_regiontseq,selist_on_nt];
                        end
                    end
                end
                CDSs.regionstr{CDSs.num,1}=regions.inforstr{regions.num,1};
                CDSs.locustag{CDSs.num,1}=CDS_locustr;
                CDSs.translation{CDSs.num,1}=CDS_translation;
                CDSs.inforstr{CDSs.num,1}=CDS_inforstr;
                CDSs.iscomplement(CDSs.num,1)=CDSiscomplement;
                CDSs.borders(CDSs.num,:)=CDS_border;
                CDSs.seqlen(CDSs.num,1)=CDS_length;
                CDSs.seposis_matrix{CDSs.num,1}=seposis_matrix;
                CDSs.seposis_sepnum(CDSs.num,1)=seposis_sepnum;
                CDSs.id_in_regiontseq{CDSs.num,1}=id_in_regiontseq;
                CDSs.seq_in_transorder{CDSs.num,1}=seq_in_transorder;
                CDSs.id_in_CDStrans{CDSs.num,1}=id_in_CDStrans;% the order of translation         
                CDSs.function{CDSs.num,1}=table_cds.gene_kind{CDS_locs_id};
                CDSs.product{CDSs.num,1}=check_field(table_cds,'product',CDS_locs_id);
                CDSs.gene_name{CDSs.num,1}=check_field(table_cds,'gene',CDS_locs_id);
                CDSs.protein_id{CDSs.num,1}=check_field(table_cds,'protein_id',CDS_locs_id);
                CDSs.region_ids(CDSs.num,1)=regions.num;
                if abs(length(id_in_regiontseq)-length(id_in_CDStrans))>4 || abs(length(seq_in_transorder)-length(id_in_CDStrans))>4
                    CDSs.translationError(CDSs.num,1)=1; % record if translation has errors
                    CDSs.translationError_length(CDSs.num,1)=abs(length(id_in_regiontseq)-length(id_in_CDStrans));
                    CDSs.translationError_length(CDSs.num,2)=abs(length(seq_in_transorder)-length(id_in_CDStrans));
                    fprintf('Translation has errors in %s\n',CDS_inforstr)% translation have errors
                else
                    CDSs.translationError(CDSs.num,1)=0;
                    CDSs.translationError_length(CDSs.num,1)=0;
                    CDSs.translationError_length(CDSs.num,2)=0;
                end
                domains_overlappercent_thisCDS=(1+min(CDS_border(2),domains_border(:,2))-max(CDS_border(1),domains_border(:,1)))./domains_ntseqlen;  % which domains overlap with this CDS
                % get the max coverage of domains on CDS
                for d=1:size(domains_maxoverlapp_cdsid)
                    if domains_overlappercent_thisCDS(d)>domains_maxoverlapp_cdsid(d,1)
                        domains_maxoverlapp_cdsid(d,1)=domains_overlappercent_thisCDS(d);
                        domains_maxoverlapp_cdsid(d,2)=CDSs.num;
                    end
                end
            end
        end
        if omains_flag
            uniCDSs_relaventtodomain=unique(domains_maxoverlapp_cdsid(:,2));
            for c=1:length(uniCDSs_relaventtodomain)
                % which domains overlap with this one;
                cid=uniCDSs_relaventtodomain(c);
                if cid>0
                    dlocs_thisc=find(domains_maxoverlapp_cdsid(:,2)==cid);
                    interdomains_border_thiscDS=zeros(length(dlocs_thisc)+1,2);
                    if CDSs.iscomplement(cid)==1 % on the reverse direction
                        [~,order]=sort(domains_border(dlocs_thisc,1),'descend');
                        interdomains_border_thiscDS(1,2)=CDSs.borders(cid,2);
                        interdomains_border_thiscDS(end,1)=CDSs.borders(cid,1);
                    else
                        [~,order]=sort(domains_border(dlocs_thisc,1),'ascend');
                        interdomains_border_thiscDS(1,1)=CDSs.borders(cid,1);
                        interdomains_border_thiscDS(end,2)=CDSs.borders(cid,2);
                    end
                    dlocs_thisc=dlocs_thisc(order);
                    CDS_trans=CDSs.translation{cid};
                    id_in_regiontseq=CDSs.id_in_regiontseq{cid};
                    seq_in_transorder=CDSs.seq_in_transorder{cid};
                    id_in_CDStrans=CDSs.id_in_CDStrans{cid};% the order of translation
                    %%%%%%%load domains and interdomains
                    % border of interdomains
                    for d=1:length(dlocs_thisc)
                        if CDSs.iscomplement(cid)==1 % reverse complement
                            interdomains_border_thiscDS(d,1)=domains_border(dlocs_thisc(d),2)+1;
                            interdomains_border_thiscDS(d+1,2)=domains_border(dlocs_thisc(d),1)-1;
                        else % forward
                            interdomains_border_thiscDS(d,2)=domains_border(dlocs_thisc(d),1)-1;
                            interdomains_border_thiscDS(d+1,1)=domains_border(dlocs_thisc(d),2)+1;
                        end
                    end
                    for d=1:(length(dlocs_thisc)+1)
                        for itd=0:1% 0: interdomain, 1: interdomain
                            semat=[];
                            if itd==1 && d~=(length(dlocs_thisc)+1)% domain
                                dloc=dlocs_thisc(d);
                                semat=domains_sematrix{dloc};
                            elseif itd==0 % interdomain
                                if d==1
                                    dtypestrbefore='Head';
                                    dtypebefore=-1;
                                else
                                    dtypestrbefore=domains_typestr{dlocs_thisc(d-1)};
                                    dtypebefore=domains_typeid(dlocs_thisc(d-1));
                                end
                                if d==(length(dlocs_thisc)+1)
                                    dtypestrafter='Tail';
                                    dtypeafter=-1;
                                else
                                    dtypestrafter=domains_typestr{dlocs_thisc(d)};
                                    dtypeafter=domains_typeid(dlocs_thisc(d));
                                end
                                semat=interdomains_border_thiscDS(d,:);
                            end
                            % obtain translation of this section
                            if ~isempty(semat)
                                % get translation of this domain
                                omain_trans_ntseq='';
                                omain_trans_aaids=[];
                                omain_trans_locs=[];
                                for s=1:size(semat,1)
                                    sp=semat(s,1);
                                    ep=semat(s,2);
                                    for r=sp:ep
                                        locinseq=find(id_in_regiontseq==r,1);
                                        if ~isempty(locinseq) && ~(length(seq_in_transorder)<locinseq) && ~(length(id_in_CDStrans)<locinseq)
                                            omain_trans_locs(end+1)=locinseq;
                                            omain_trans_ntseq(end+1)=seq_in_transorder(locinseq);
                                            omain_trans_aaids(end+1)=id_in_CDStrans(locinseq);
                                        end
                                    end
                                end
                                [~,order]=sort(omain_trans_locs);
                                omain_trans_ntseq=omain_trans_ntseq(order);
                                omain_trans_aaids=omain_trans_aaids(order);
                                omain_trans_aaseq=CDS_trans(unique(omain_trans_aaids));
                                % add it to domain or interdomain
                                omains.num=omains.num+1;
                                if itd==0 % interdomain
                                    omains.iscomplement(omains.num,1)=CDSs.iscomplement(cid);
                                    omains.subtypestr{omains.num,1}='nan';
                                    omains.locustag{omains.num,1}= CDSs.locustag{cid};
                                    omains.translation{omains.num,1}='nan';
                                    omains.typeid_mat(omains.num,:)=[dtypebefore,dtypeafter];
                                    omains.typestr{omains.num,1}=['Between_',dtypestrbefore,'_',dtypestrafter];
                                    omains.specificity{omains.num,1}=[];
                                else% domain
                                    omains.iscomplement(omains.num,1)= domains_iscomplement(dloc);
                                    omains.subtypestr{omains.num,1}=domains_subtypestr{dloc};
                                    omains.locustag{omains.num,1}= domains_locustag{dloc};
                                    omains.translation{omains.num,1}=domains_translation{dloc}; % translation annotated by antismash which is always short in A and C domain
                                    omains.typeid_mat(omains.num,:)=[domains_typeid(dloc),0];
                                    omains.typestr{omains.num,1}=domains_typestr{dloc};
                                    omains.specificity{omains.num,1}=domains_specificity_mat{dloc};
                                end
                                omains.CDSstr{omains.num,1}=CDSs.inforstr{cid};
                                omains.seq_ntaa{omains.num,1}=omain_trans_ntseq;
                                omains.seq_ntaa{omains.num,2}=omain_trans_aaseq;
                                omains.isdomain(omains.num,1)=itd;
                                omains.borders(omains.num,:)=[semat(1),semat(end)];
                                omains.inforstr{omains.num,1}=[CDSs.inforstr{cid},'_',num2str(d),'_',omains.typestr{omains.num,1}];
                                omains.CDS_ids(omains.num,1)=cid; % 指的是domain的CDS id, interdomain可能横跨多个CDS
                                omains.region_ids(omains.num,1)=regions.num;
                            end
                        end
                    end
                end
            end
        end
    end
end
% toc
if regions_flag
    rawset.regions=regions;
end
if cand_clusters_flag
    rawset.cand_clusters=cand_clusters;
end
if CDSs_flag
    rawset.CDSs=CDSs;
end
if omains_flag
    rawset.omains=omains;
end
end
%% 
function loc_str=strremovemulti_split(loc_str,strremove_list,delimiter)
for i = 1:length(strremove_list)
    loc_str=strrep(loc_str,strremove_list{i},'');
end
loc_str=split(loc_str,delimiter);
end

function [CDS_locustr,CDS_translation,CDSiscomplement,CDS_border,seposis_sepnum,seposis_matrix,CDS_length] = infor_generate(table_cds,CDS_locs_id)
if ismember('locus_tag',table_cds.Properties.VariableNames)%table_cds是table不是struct，不能用isfield来判断是否有有特定字段
    CDS_locustr=table_cds.locus_tag{CDS_locs_id};
elseif ismember('gene',table_cds.Properties.VariableNames)
    CDS_locustr=table_cds.gene{CDS_locs_id}; % 有些基因组注释把locus_tag放在gene name里了
elseif ismember('protein_id',table_cds.Properties.VariableNames)
    CDS_locustr=table_cds.protein_id{CDS_locs_id}; % 实在不行用protein_id
else
    CDS_locustr=[];% 都没有就留空
end
CDS_translation=table_cds.translation{CDS_locs_id};
if strcmp(table_cds.Direction{CDS_locs_id},'forward')
    CDSiscomplement=0;
else
    CDSiscomplement=1;
end
CDS_border=[table_cds.Minimum(CDS_locs_id),table_cds.Maximum(CDS_locs_id)];
loc_location=strremovemulti_split(table_cds.Location{CDS_locs_id},{'complement','join','(',')','<','>'},',');
seposis_sepnum=length(loc_location);
seposis_matrix=zeros(seposis_sepnum,2);
for i = 1:seposis_sepnum
    seposis_matrix(i,:)=str2double(split(loc_location{i},'..'));
end
CDS_length=table_cds.Length(CDS_locs_id);
end

function intersection_length = cal_inter_len(range1,range2)
% 计算交集范围
intersection_range = [max(range1(1), range2(1)), min(range1(2), range2(2))];

% 计算交集长度
intersection_length = max(0, intersection_range(2) - intersection_range(1) + 1);
end

function max_core_overlap_ratio=cal_max_overlap_ratio(Cand_cluster,proto_core,cmp_proto_core_id,cmp_Cand_cluster_id,protocluster)
% input没protocluster的时候，算的是core;有protocluster的时候，算的是neighbourhood
cmp_proto_core=str2double(Cand_cluster.protoclusters{cmp_proto_core_id});
if nargin == 4
    loc_core_overlap_ratio_list=zeros(length(cmp_proto_core),1);
elseif nargin == 5
    loc_core_overlap_ratio_list=zeros(length(cmp_proto_core),2);
end
for core_id = 1:length(cmp_proto_core)
    loc_index_core=proto_core.protocluster_number==cmp_proto_core(core_id);
    if nargin == 4
        loc_core_overlap_ratio_list(core_id)=cal_inter_len([proto_core.Minimum(loc_index_core),proto_core.Maximum(loc_index_core)],[Cand_cluster.Minimum(cmp_Cand_cluster_id),Cand_cluster.Maximum(cmp_Cand_cluster_id)])/proto_core.Length(loc_index_core);
    elseif nargin == 5
        loc_index_cluster=protocluster.protocluster_number==cmp_proto_core(core_id);
        loc_core_overlap_ratio_list(core_id,1)=cal_inter_len([protocluster.Minimum(loc_index_cluster),proto_core.Minimum(loc_index_core)],[Cand_cluster.Minimum(cmp_Cand_cluster_id),Cand_cluster.Maximum(cmp_Cand_cluster_id)])/(proto_core.Minimum(loc_index_core)-protocluster.Minimum(loc_index_cluster)+1);
        loc_core_overlap_ratio_list(core_id,2)=cal_inter_len([proto_core.Maximum(loc_index_core),protocluster.Maximum(loc_index_cluster)],[Cand_cluster.Minimum(cmp_Cand_cluster_id),Cand_cluster.Maximum(cmp_Cand_cluster_id)])/(protocluster.Maximum(loc_index_cluster)-proto_core.Maximum(loc_index_core)+1);
    end
end
max_core_overlap_ratio=max(loc_core_overlap_ratio_list,[],'all');
end

function proto_core=extract_infor(proto_core,table_other,Category_str)
loc_inex=ismember(table_other.Category,Category_str);
proto_core.protocluster_number=str2double(table_other.protocluster_number(loc_inex));% id
proto_core.Minimum=table_other.Minimum(loc_inex);
proto_core.Maximum=table_other.Maximum(loc_inex);
end

function [cluster_length,cds_length] = calculate_BGC_length(table_cds,CDS_locs_id_list,cds_type)
if strcmp(cds_type,'biosynthetic')
    loc_gene_kind_index=ismember(table_cds.gene_kind(CDS_locs_id_list),'biosynthetic');
elseif strcmp(cds_type,'biosynthetic+')
    loc_gene_kind_index=ismember(table_cds.gene_kind(CDS_locs_id_list),'biosynthetic')|ismember(table_cds.gene_kind(CDS_locs_id_list),'biosynthetic-additional');
elseif strcmp(cds_type,'functional')
    loc_gene_kind_index=~ismember(table_cds.gene_kind(CDS_locs_id_list),'other');
elseif strcmp(cds_type,'all')
    loc_gene_kind_index=1:length(CDS_locs_id_list);
end
cluster_length=max(table_cds.Maximum(CDS_locs_id_list(loc_gene_kind_index)))-min(table_cds.Minimum(CDS_locs_id_list(loc_gene_kind_index)))+1;
cds_length=sum(table_cds.Length(CDS_locs_id_list(loc_gene_kind_index)));
end

function result = check_field(table_cds,field_str,CDS_locs_id)
if ismember(field_str,table_cds.Properties.VariableNames)
    result=table_cds.(field_str){CDS_locs_id};
else
    result=[];
end
end