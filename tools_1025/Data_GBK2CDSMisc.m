function [feature_cds_trim,feature_misc_trim,feature_gene_trim,feature_cds_trash,tempgbk] = Data_GBK2CDSMisc(gbk)
% extracts gene information about CDS and other genes (tRNA,rRNA etc.) from
% gbk file and saves their info to structs separately
% INPUT ------------------------------------------------------------
% gbkadrs: a gbk file or file address of the gbk file
% OUTPUT ----------------------------------------------------------
% feature_cds_trim: struct storing information about CDS in the genome
% feature_misc_trim: struct storing information about other genes in the
% genome, including tRNA, rRNA etc.
% feature_gene_trim: struct storing information about gene in the genome
% feature_cds_trash: struct storing information about CDS containing unwanted tags
% tempgbk: result of genbankread_HRL(gbk)
%% Initialization
% Loading data
if ischar(gbk)
    tempgbk = genbankread_HRL(gbk); % to handle no feature gbk
elseif isstruct(gbk)
    tempgbk = gbk;
    clear gbk
end
gbk_index=[];
tempfeature_raw=[];% 解决了tempgbk有不止一条序列时的问题。之前是tempfeature = tempgbk.Features;如果tempgbk是n*1的struct，则只能读到最后一个的内容。
for i = 1:length(tempgbk)
    if ~isempty(tempgbk(i).Features) % to handle no feature gbk
        tempfeature_raw=[tempfeature_raw;cellstr(tempgbk(i).Features)];
        gbk_index=[gbk_index;repmat(i,size(tempgbk(i).Features,1),1)]; % record each feature from which gbk file
    end
end
Col_n_list=cellfun(@numel, tempfeature_raw);
maxCols=max(Col_n_list);
tempfeature_raw = cellfun(@(x) [x, repmat(' ', 1, maxCols - numel(x))], tempfeature_raw, 'UniformOutput', false);
tempfeature = char(tempfeature_raw);
if ~isempty(tempfeature)
    %% For GenBank Contigs
    if ~isfield(tempgbk,'CDS')
        for i = 1:length(tempgbk)
            tempgbk(i).CDS = {};
        end
    end
    %% Features to Cat&Comment struct
    % Split field name and field value by its format in the feature text
    feature = struct('Category',{},'Comment',{});
    num = 0;
    gbk_index_list=[]; % record each feature from which gbk file
    % location infos
    cmt = tempfeature(1,17:end);
    cat = strtrim(tempfeature(1,1:16));
    % Check the data from the second line (because first line is the location info)
    for lineidx = 2:size(tempfeature,1)
    % there are 16 space in lines_cmt if it isn't a head
        tempcat = tempfeature(lineidx,1:16);
        tempcmt = tempfeature(lineidx,17:end);
        % A line following the precedent
        if startsWith(tempcat,' ')
            cmt = [cmt;tempcmt];
        % A line starts a new entry
        else
            gbk_index_list=[gbk_index_list;gbk_index(lineidx-1)];
            num = num+1;
            feature(num).Category = cat;
            feature(num).Comment = cmt;
            
            cat = strtrim(tempcat);
            cmt = tempcmt;
        end
    end
    gbk_index_list=[gbk_index_list;gbk_index(lineidx)];
    num = num+1;
    feature(num).Category = cat;
    feature(num).Comment = cmt;
    feature_nogene = feature;
    %% Convert comment to struct
    feature_trim = feature_nogene;
    feature_cds = struct('Minimum',{},'Maximum',{},'Length',{},'Direction',{},'Comment',{});
    feature_gene = struct('Minimum',{},'Maximum',{},'Length',{},'Direction',{},'Comment',{});
    feature_misc = struct('Category',{},'Minimum',{},'Maximum',{},'Length',{},'Direction',{},'Comment',{});
    cdsnum = 0;
    genenum = 0;
    miscnum = 0;
    for ftidx = 1:length(feature_nogene)
        cat = feature_nogene(ftidx).Category;
        cmt = feature_nogene(ftidx).Comment;
        x=tempgbk(gbk_index_list(ftidx));          
        switch cat
            case 'gene'
                genenum = genenum + 1;
                feature_gene = Data_AddFieldInfo(cmt,feature_gene,genenum,x.Accession);
                feature_trim = Data_AddFieldInfo(cmt,feature_trim,ftidx,x.Accession);
            case 'CDS'
                cdsnum = cdsnum + 1;
                feature_cds = Data_AddFieldInfo(cmt,feature_cds,cdsnum,x.Accession);
                feature_trim  = Data_AddFieldInfo(cmt,feature_trim,ftidx,x.Accession);
            otherwise
                miscnum = miscnum + 1;
                feature_misc = Data_AddFieldInfo(cmt,feature_misc,miscnum,x.Accession);
                feature_trim = Data_AddFieldInfo(cmt,feature_trim,ftidx,x.Accession); 
                feature_misc(miscnum).Category = cat;  
        end
    end
    feature_cds_trim = feature_cds;
    feature_cds_trash = struct('Minimum',[]);
    feature_misc_trim = feature_misc;
    feature_gene_trim = feature_gene;
else
    feature_cds_trim = struct();
    feature_misc_trim = struct();
    feature_gene_trim = struct();
    feature_cds_trash = struct();
end
end   

function feature = Data_AddFieldInfo(comment,feature,idx,Accession)
% feature = Data_AddFieldInfo(comment,feature,idx)
% Categorizes info in 'comment' and adds them to different fields of a
% struct array 'feature' separately. If that field does not exist, a new
% field will be created.
% INPUT --------------------------------------------
% comment: character array, the first line is the location info
%          new entries begin with '/', field names and values are separated
%          by '='
% feature: the feature struct array to be modified (i.e. to which new
%          info/fields to be added)
% idx: the index for the target struct in 'feature'
% EXAMPLE ------------------------------------------
% feature = Data_AddFieldInfo(comment,feature,3);
% NOTE --------------------------------------------
% Based on the genbank Comment format, the name of a non-empty field will
% be followed by a equal sign, behind which is the field content.
% Character content is enclosed by double quotation marks.
% e.g. numeric field: '/codon_start=1                               '
%      string field: '/protein_id="BAO74231.1"                      '
%      multiline field: '/translation="MEITAQSVWSNCLAFIKDNIQSQAYKTWFEPIDA'
%                         'VKLSGNALSIQVPSKFFYEWLEEHYVKILKVSLTKELGEDAKLVYVIK'
%                         'MENTYGNKQPFTEKIPSSNR"                           '
%      tag field: '/pseudo                                        '

%% Add values to necessary fields
feature(idx).Comment = comment;
start = find(comment(:,1) == '/',1);
if isempty(start)
    location_last_line=size(comment,1);
    new_comment = '';
else
    location_last_line=start-1;
    new_comment = comment(start:end,:);
end
loc_location=[];% location could be very long, in the more than one lines
for i = 1:location_last_line
    loc_location=[loc_location,comment(i,:)];
end
comment=new_comment;
[Max,Min,Len,Dir] = Data_LocateCDS(loc_location);
feature(idx).Location = loc_location;

feature(idx).Minimum = Min;
feature(idx).Maximum = Max;
feature(idx).Length = Len;
feature(idx).Direction = Dir;


%% Iteration to add optional fields

partialopt = 0;  
anticodonopt = 0;

for lineidx = 1:size(comment,1)
    newline = comment(lineidx,:);
    newline = strtrim(newline);

    % All the precedent fields are completed
    if ~partialopt
        equpin = find(newline == '=');

        % for fields with content, e.g. /note
        if ~isempty(equpin)

            switch sum(newline == '"')
                % for numeric content
                case 0
                    field = newline(2:equpin-1);

                    % for anticodon field
                    if strcmp(field,'anticodon')
                        anticodonopt = 1;
                        partialopt = 1;
                        val = newline(equpin+1:end);
                        if lineidx==size(comment,1)||contains(strtrim(comment(lineidx+1,:)),'=') % check anticodon if is completed
                            feature = Feature_update(feature,idx,field,val);
                            anticodonopt = 0;
                            partialopt = 0;
                        end
                    else
                        val = str2double(newline(equpin+1:end));
                        feature = Feature_update(feature,idx,field,val); 
                    end
                % for complete character content    
                case 2
                    field = newline(2:equpin-1);
                    val = newline(equpin+2:end-1);
                    feature = Feature_update(feature,idx,field,val);     
                % for partial character content     
                case 1
                    field = newline(2:equpin-1);
                    val = newline(equpin+2:end);
                    partialopt = 1;
            end

        % for tag fields, e.g. /pseudo
        else
            field = newline(2:end);
            val = 1;

            feature = Feature_update(feature,idx,field,val);
        end

    % the current field is waiting to be completed
    else
        % the current field is completed
        if strcmp(field,'translation') % no space in translation
            loc_delimiter='';
        else
            loc_delimiter=' ';
        end
        if contains(newline,'"')
            val = [val loc_delimiter newline(1:end-1)];
            feature = Feature_update(feature,idx,field,val);
            partialopt = 0;
        % complete anticodon field    
        elseif anticodonopt
            val = [val loc_delimiter newline];
            feature = Feature_update(feature,idx,field,val);
            anticodonopt = 0;
            partialopt = 0;
        % the current field is still unfinished
        else
            val = [val loc_delimiter newline];
        end  
    end
end

if true
    feature(idx).Accession = Accession;
end
end

function feature = Feature_update(feature,idx,field,val)
% for handle the feature which has more than one items
% note: when the type of items is "char", if they only contains one items, it
% will be "char", but it will be "cell" if it has more than one items
    if isfield(feature(idx),field)&&~isempty(feature(idx).(field))
        if ischar(val)
            feature(idx).(field) = [feature(idx).(field);{val}];
        else
            feature(idx).(field) = [feature(idx).(field);val];
        end
    else
        feature(idx).(field) = val;
    end
end

function [Max,Min,Len,Dir] = Data_LocateCDS(cdsinfo)
% extracts information about a cds, including the maximum, minimum,
% direction, and length. 
% INPUT ----------------------------------------------
% cdsinfo: string array, NCBI format, e.g. 'complement(43954..45036)  '

if contains(cdsinfo, 'complement')
    Dir = 'reverse';
    cdsinfo = cdsinfo(find(cdsinfo == '(')+1 : find(cdsinfo == ')')-1);
else
    Dir = 'forward';
end

match = regexp(cdsinfo,'[0-9]*','Match');
match = cellfun(@str2double,match);
Min = min(match);
Max = max(match);

Len = Max - Min + 1;

end