function data=genbankread_HRL(gbtext,NVPArgs)
%GENBANKREAD reads GenBank format data files.
%
%   DATA = GENBANKREAD(FILE) reads in a GenBank formatted sequence from
%   FILE and creates a structure DATA containing fields corresponding to
%   the GenBank keywords. If the file contains information about multiple
%   sequences, then the information is stored in an array of
%   structures.
%
%   FILE is a character vector or string specifying a file name or a URL.
%   It can also be a character array or string vector that contains the
%   text of a GenBank format file.
%
%   GENBANKREAD(URL, 'Timeout', T) sets the connection timeout to read
%   the GenBank file located at URL to T, in seconds. Default is 5.
%
%   Based on version 179.0 of GenBank
%
%   Examples:
%
%       % Download a GenBank file to your local drive.
%       getgenbank('M10051', 'TOFILE', 'HGENBANKM10051.GBK')
%
%       % Then bring it into a MATLAB sequence.
%       data = genbankread('HGENBANKM10051.GBK')
%
%   See also EMBLREAD, FASTAREAD, GENPEPTREAD, GETGENBANK, SCFREAD, 
%   SEQVIEWER.

% Copyright 2002-2020 The MathWorks, Inc.
%
% Modified by Ruolin He in 2023/12/30 to handle no feature gbk
% if ~matchstart(gbtext(ln,:),'ORIGIN') && ~matchstart(gbtext(ln,:),'CONTIG') % handle no feature gbk
% if ~isempty(data(record_count).Features)% handle no feature gbk
%
% Modified in function extract_feature by Ruolin He in 2023/12/30 to avoid wrong match
% such as "CDS_motif"
% theFeature=[theFeature,' ']; % add a blank to avoid wrong match
%
% Modified in function extract_feature by Ruolin He in 2023/12/30 to solve
% buffer error (bufsize = 4095) when translation sequence is large
% theStruct(featloop).translation = strcat(fullstr{:})
%
% Modified in function extract_feature by Ruolin He in 2023/12/30 to solve
% if the first line in Comment is blank, comment will be empty
% data(record_count).Comment=strvcat(data(record_count).Comment, strtrim(gbtext(ln,1:t{1}(2))));

arguments
    gbtext
    NVPArgs.Timeout(1,1) double {mustBeNumeric, mustBePositive} = 5;
    NVPArgs.PreambleText(1,1) logical  = false;
    NVPArgs.AllFeatures(1,1) logical = false;
end

gbtext = convertStringsToChars(gbtext);

getPreambleText = NVPArgs.PreambleText;
allFeatures = NVPArgs.AllFeatures;

if ~ischar(gbtext) && ~iscellstr(gbtext)
    error(message('bioinfo:genbankread:InvalidInput'));
end

if iscellstr(gbtext)
    % do not mess with it, just put it to char and try to decipher in the try-catch trap
    gbtext=char(gbtext);
else %it is char, lets check if it has an url or a file before try to decipher it
    if size(gbtext,1)==1 && ~isempty(strfind(gbtext(1:min(10,end)), '://'))
        % must be a URL
        try
            options = weboptions('ContentType','text', 'Timeout', NVPArgs.Timeout);
            gbtext = webread(gbtext,options);
        catch allExceptions
            error(message('bioinfo:genbankread:CannotReadURL', gbtext));
        end
        % clean up any &amp s
        gbtext=strrep(gbtext,'&amp;','&');
        % make each line a separate row in a string array
        gbtext = textscan(gbtext,'%s','delimiter','\n','whitespace','');
        gbtext = char(gbtext{1});
    elseif size(gbtext,1)==1
        if (exist(gbtext,'file') || exist(fullfile(pwd,gbtext),'file'))
            fid = fopen(gbtext,'r');
            gbtext = textscan(fid,'%s','delimiter','\n','whitespace','');
            gbtext = char(gbtext{1});
            fclose(fid);   
        else
            gbtext = textscan(gbtext,'%s','delimiter','\n','whitespace','');
            gbtext = char(gbtext{1});
        end
    end
end

% If the input is a string of GenBank data then words LOCUS and DEFINITION must be present
if size(gbtext,1)==1 || isempty(strfind(gbtext(1,:),'LOCUS')) || isempty(strfind(gbtext(2,:),'DEFINITION'));
    error(message('bioinfo:genbankread:NonMinimumRequiredFields'))
end

%line number
ln = 1;

%multiple records possible in one record
record_count=1;

numLines = size(gbtext,1);
try
    while 1,

        %LOCUS - Mandatory keyword/exactly one record.
        loc_locus=split(gbtext(ln,:));
        data(record_count).LocusName = loc_locus{2};  %#ok<*AGROW>
        data(record_count).LocusSequenceLength =loc_locus{3};
        data(record_count).LocusNumberofStrands = strtrim(gbtext(ln,45:47));
        data(record_count).LocusTopology = strtrim(gbtext(ln,56:63));
        data(record_count).LocusMoleculeType = strtrim(gbtext(ln,48:53));
        data(record_count).LocusGenBankDivision = strtrim(gbtext(ln,65:67));
        data(record_count).LocusModificationDate = strtrim(gbtext(ln,69:79));

        ln=ln+1;

        %DEFINITION - Mandatory keyword/one or more records.
        [~,~,t] = regexp(gbtext(ln,:),'DEFINITION\s+(\w|\W)+'); 
        data(record_count).Definition = strtrim(gbtext(ln,t{1}(1):t{1}(2)));

        ln=ln+1;

        while ~matchstart(gbtext(ln,:),'ACCESSION')
            data(record_count).Definition = [data(record_count).Definition,' ', strtrim(gbtext(ln,:))];
            ln = ln+1;
        end

        %ACCESSION - Mandatory keyword/one or more records.
        [~,~,t] = regexp(gbtext(ln,:),'ACCESSION\s+(\w|\W)+'); 
        data(record_count).Accession = strtrim(gbtext(ln,t{1}(1):t{1}(2)));
        ln=ln+1;

        while ~matchstart(gbtext(ln,:),'VERSION')
            data(record_count).Accession=[data(record_count).Accession ' ' strtrim(gbtext(ln,:))];
            ln=ln+1;
        end

        %VERSION - Mandatory keyword/exactly one record.
        [~,~,t] = regexp(gbtext(ln,:),'VERSION\s+(\w|\W)+'); 
        rmdr = gbtext(ln,t{1}(1):t{1}(2));
        [data(record_count).Version, rmdr] = strtok(rmdr,'GI:'); %#ok<STTOK>
        data(record_count).Version = strtrim(data(record_count).Version);

        %GI - Mandatory (part of version)
        data(record_count).GI = deblank(rmdr(4:end));

        ln=ln+1;

        % NID - Obsolete since 12/1999
        
        %PROJECT - Introduced by NCBI in Feb 2006, treated as optional.
        %          Obsoleted in Release 171.0 April 2009
        data(record_count).Project=[];
        [s,~,t] = regexp(gbtext(ln,:),'PROJECT\s+(\w|\W)+'); 
        if ~isempty(s)
            data(record_count).Project=deblank(gbtext(ln,t{1}(1):t{1}(2)));
            ln=ln+1;
        end
        
        %DBLINK - New field introduced by NCBI in Feb 2009, treated as
        %         Optional keyword/one or more records.
        data(record_count).DBLink = [];
        [s,~,t] = regexp(gbtext(ln,:),'DBLINK\s+(\w|\W)+'); 
        if ~isempty(s)
            data(record_count).DBLink=deblank(gbtext(ln,t{1}(1):t{1}(2)));
            ln=ln+1;
        end
        
        while ~matchstart(gbtext(ln,:),'KEYWORDS')
            data(record_count).DBLink=[data(record_count).DBLink ' ' strtrim(gbtext(ln,:))];
            ln=ln+1;            
        end

        %KEYWORDS - Mandatory keyword in all annotated entries/one or more
        %           records. 
        data(record_count).Keywords=[];
        [s,~,t] = regexp(gbtext(ln,:),'KEYWORDS\s+(\w|\W)+'); 
        if ~isempty(s)
            data(record_count).Keywords=deblank(gbtext(ln,t{1}(1):t{1}(2)));
            ln=ln+1;
        end
        while ~isempty(s) && ~matchstart(gbtext(ln,:),'SEGMENT') && ~matchstart(gbtext(ln,:),'SOURCE') 
            data(record_count).Keywords=strvcat(data(record_count).Keywords, deblank(gbtext(ln,:))); %#ok<*VCAT>
            ln=ln+1;
        end
        if all(~isletter(data(record_count).Keywords)),
            data(record_count).Keywords = [];
        end


        %SEGMENT -  Optional keyword (only in segmented entries)/exactly
        %           one record. 

        data(record_count).Segment=[];
        [s,~,t] = regexp(gbtext(ln,:),'SEGMENT\s+(\w|\W)+'); 
        if ~isempty(s)
            data(record_count).Segment=gbtext(ln,t{1}(1):t{1}(2));
            ln=ln+1;
        end


        %SOURCE - Mandatory keyword in all annotated entries/one or
        %         more records/includes one subkeyword.
        [~,~,t] = regexp(gbtext(ln,:),'SOURCE\s+(\w|\W)+'); 
        data(record_count).Source = deblank(gbtext(ln,t{1}(1):t{1}(2)));
        ln=ln+1;
        while ~matchstart(gbtext(ln,:),'ORGANISM') && ~matchstart(gbtext(ln,:),'FEATURES') && ~matchstart(gbtext(ln,:),'COMMENT') && ~matchstart(gbtext(ln,:),'BASE COUNT')
            data(record_count).Source = [data(record_count).Source ' ' deblank(gbtext(ln,t{1}(1):t{1}(2)))];
            ln=ln+1;
        end


        %ORGANISM - Mandatory for all annotated records.
        data(record_count).SourceOrganism = [];
        [s,~,t] = regexp(gbtext(ln,:),'ORGANISM\s+(\w|\W)+'); 
        if ~isempty(s)
            data(record_count).SourceOrganism = strtrim(gbtext(ln,t{1}(1):t{1}(2)));
            ln=ln+1;
            while ~matchstart(gbtext(ln,:),'REFERENCE') && ~matchstart(gbtext(ln,:),'COMMENT')
                data(record_count).SourceOrganism = strvcat(data(record_count).SourceOrganism, strtrim(gbtext(ln,:))); 
                ln=ln+1;
            end
        end

        %REFERENCE
        [data,gbtext,ln] = referenceparse(data,gbtext,ln,record_count);

        %COMMENT - Optional
        data(record_count).Comment = [];
        [s,~,t] = regexp(gbtext(ln,:),'COMMENT\s+(\w|\W)+'); 
        if ~isempty(s)
            data(record_count).Comment = strtrim(gbtext(ln,t{1}(1):t{1}(2)));
            ln=ln+1;
            while ~matchstart(gbtext(ln,:),'FEATURES') && ~matchstart(gbtext(ln,:),'BASE COUNT')...
                    && ~matchstart(gbtext(ln,:),'ORIGIN')
%                 data(record_count).Comment=strvcat(data(record_count).Comment, strtrim(gbtext(ln,t{1}(1):t{1}(2)))); 
                data(record_count).Comment=strvcat(data(record_count).Comment, strtrim(gbtext(ln,1:t{1}(2)))); % solve if the first line in Comment is blank
                ln=ln+1;
            end
        end

        %FEATURES - Optional
        data(record_count).Features = [];
        lnFeatures = inf;
        if matchstart(gbtext(ln,:),'FEATURES')
            lnFeatures = ln; % save this position to get preamble text later
            feats = cell(numLines-ln,1);
            ln=ln+1;
            if ~matchstart(gbtext(ln,:),'ORIGIN') && ~matchstart(gbtext(ln,:),'CONTIG') % handle no feature gbk
                featCount = 1;
                feats{featCount} = gbtext(ln,1:end);
                ln=ln+1;
                while ~matchstart(gbtext(ln,:),'ORIGIN') && ~matchstart(gbtext(ln,:),'CONTIG')
                    featCount = featCount+1;
                    feats{featCount} = gbtext(ln,1:end);
                    ln=ln+1;
                end
                data(record_count).Features = strtrim(strvcat(feats(1:featCount))); 
            end
            % Extract information for the features -- CDS, gene, mRNA,
            % misc_feature and repeat_unit
            if ~isempty(data(record_count).Features)% handle no feature gbk
                try
                    data(record_count).CDS = extract_feature(data(record_count).Features,'CDS');
                catch allExceptions %#ok<NASGU>
                    if isfield(data(record_count),'CDS')
                        data(record_count).CDS = [];
                    end
                    warning(message('bioinfo:genbankread:BADCDS'));
                end
            end
            if allFeatures
                try
                    data(record_count).mRNA = extract_feature(data(record_count).Features,'mRNA');
                catch allExceptions %#ok<NASGU>
                    if isfield(data(record_count),'mRNA')
                        data(record_count).mRNA = [];
                    end
                    warning(message('bioinfo:genbankread:BADmRNA'));
                end
                try
                    data(record_count).gene = extract_feature(data(record_count).Features,'gene');
                catch alLExceptions %#ok<NASGU>
                    if isfield(data(record_count),'gene')
                        data(record_count).gene = [];
                    end
                    warning(message('bioinfo:genbankread:BADgene'));
                end
                try
                    data(record_count).misc_feature = extract_feature(data(record_count).Features,'misc_feature');
                catch allExceptions %#ok<NASGU>
                    if isfield(data(record_count),'misc_feature')
                        data(record_count).misc_feature = [];
                    end
                    warning(message('bioinfo:genbankread:BADmisc_feature'));
                end
                try
                    data(record_count).repeat_unit = extract_feature(data(record_count).Features,'repeat_unit');
                catch allExceptions %#ok<NASGU>
                    if isfield(data(record_count),'repeat_unit')
                        data(record_count).repeat_unit = [];
                    end
                    warning(message('bioinfo:genbankread:BADrepeat_unit'));
                end
            end
        end

        %BASECOUNT - obsolete, removed October 2003
        
        %ORIGIN - Mandatory
        if matchstart(gbtext(ln,:),'ORIGIN')
            lnOrigin = ln; % save this position to get preamble text later
            ln=ln+1;
        end

		if matchstart(gbtext(ln,:),'CONTIG')
			% CONTIG
			data(record_count).Contig = '';
            startln = ln;
            % the sequence will go up to the start of the next possible record
            ln = find(gbtext(ln:end,1) == '/' | gbtext(ln:end,1) == 'L',1) + ln -1;
            if isempty(ln)
                ln = size(gbtext,1);
            end
            endln = ln-1;
			contig = strtrim(regexprep(cellstr(gbtext(startln:endln,:)),'^CONTIG',''));
            data(record_count).Contig = [contig{:}];
			warning(message('bioinfo:genbankread:contigGenBankData'));
		else	
            % SEQUENCE
            data(record_count).Sequence = '';
            startln = ln;
            % the sequence will go up to the start of the next possible record
            ln = find(gbtext(ln:end,1) == '/' | gbtext(ln:end,1) == 'L',1) + ln -1;
            if isempty(ln)
                ln = size(gbtext,1);
            end
            endln = ln-1;
            seq = gbtext(startln:endln,:)';
            seq = seq(:)';
            seq(~isletter(seq)) = '';
            data(record_count).Sequence = seq;
		end

        % PREAMBLETEXT
        if getPreambleText
            data(record_count).PreambleText = gbtext(1:min(lnOrigin,lnFeatures)-1,:);
        end

        if ln < numLines && matchstart(gbtext(ln,:),'//')
            while ln<numLines
                ln=ln+1;
                % another record ?
                if matchstart(gbtext(ln,:),'LOCUS')
                    record_count = record_count+1;
                    break
                end
            end
        end
        if ln == numLines
            return
        end
    end
catch le
    if strcmpi(le.identifier,'matlab:nomem')
        clear data
        rethrow(le)
    else
        warning(message('bioinfo:genbankread:incompleteGenBankData'));
    end
end

%-------------------------------------------------------------------------%
function theStruct = extract_feature(theText,theFeature)
% Extract the feature information
% typically CDS, gene, mRNA
theCellstr = cellstr(theText);
theFeature=[theFeature,' ']; % add a blank to avoid wrong match
% look for the tag at the start of a line
feat = strmatch(theFeature,strtrim(theCellstr)); 

featCount = numel(feat);

if featCount == 0
    theStruct = [];
    return;
end

% As the items in the CDS fields is unknown, we rely on indentation to tell
% when the CDS field ends. This is not particularly robust.
indent = find(theCellstr{feat(1)} == theFeature(1),1)-1;
if indent>0
    theCellstr = cellstr(theText(:,indent+1:end));
end
%look for lines with first level features.
featureLines = find(cellfun('isempty',regexp(theCellstr,'^\s')));
featureLines(end+1) = numel(theCellstr)+1;

% Sometimes 'CDS' show up in the /note as a single line
theFeat = intersect(featureLines,feat);
featCount= numel(theFeat);
% create empty struct
theStruct(featCount).location = '';
theStruct(featCount).gene = '';
theStruct(featCount).product = '';
theStruct(featCount).codon_start = [];
theStruct(featCount).indices = [];
theStruct(featCount).protein_id = '';
theStruct(featCount).db_xref = '';
theStruct(featCount).note = '';
theStruct(featCount).translation = '';
theStruct(featCount).text = '';

% loop over all of the CDS
for featloop = 1:featCount
    featurePos = find(featureLines == theFeat(featloop));

    endLine = featureLines(featurePos+1)-1;
    textChunk = strtrim(theCellstr(theFeat(featloop):endLine));

    numLines = numel(textChunk);
    theStruct(featloop).text = char(textChunk);

    textChunk{1} = strtrim(strrep(textChunk{1},theFeature,''));
    [featLocation, startLoop] = getFullText(textChunk, 1);
    theStruct(featloop).location = char(strread([featLocation{:}], '%s'));
    theStruct(featloop).indices = featurelocation(theStruct(featloop).location);

    strLoop = startLoop;
    while(strLoop < numLines)
        % if featloop == 618&&strLoop==23
        %     strLoop
        % end
        [fullstr, endpos] = getFullText(textChunk, strLoop);
        lines = size(fullstr, 1);
        [token,rest] = strtok(fullstr{1},'='); 
        rest= strrep(rest,'"','');
        fullstr{1}= rest(2:end);
        if(lines > 1)
            fullstr{lines} = strrep(fullstr{lines}, '"','');
        end
        switch token(2:end)
            case 'gene'
                theStruct(featloop).gene = char(fullstr{:});
            case 'product'
                theStruct(featloop).product =  char(fullstr{:});
            case 'codon_start'
                theStruct(featloop).codon_start =  char(fullstr{:});
            case 'protein_id'
                theStruct(featloop).protein_id = char(fullstr{:});
            case 'db_xref'
                theStruct(featloop).db_xref = char(fullstr{:});
            case 'note'
                theStruct(featloop).note = char(fullstr{:});
            case 'translation'
                % theStruct(featloop).translation = char(strread([fullstr{:}],'%s')); 可能会发生缓冲区错误。error: 发生缓冲区溢出(bufsize = 4095)
                theStruct(featloop).translation = strcat(fullstr{:});
            otherwise % There may be other fields...
                % disp(sprintf('Unknown field: %s',token));
        end
        strLoop = endpos;
    end
end
%-------------------------------------------------------------------------%
function [fullText, endPos] = getFullText(textChunk, i)
% next qualifier (if exists) starts with '/'
nextKey = find(strncmp('/',textChunk(i+1:end),1));
if isempty(nextKey)
    endPos = numel(textChunk)+1;
else
    endPos = i + nextKey(1);
end
fullText = textChunk(i:endPos-1);
function tf = matchstart(string,pattern)
%MATCHSTART matches start of string with pattern, ignoring spaces

% Copyright 2003-2004 The MathWorks, Inc.


tf = ~isempty(regexp(string,['^(\s)*?',pattern],'once'));

function [data,gptext,ln] = referenceparse(data,gptext,ln,record_count)
%REFERENCEPARSE parses the reference entries of GenBank files

% Copyright 2003-2006 The MathWorks, Inc.


ref_count=1;
data(record_count).Reference = [];
while ~matchstart(gptext(ln,:),'FEATURES') && ~matchstart(gptext(ln,:),'COMMENT') && ~matchstart(gptext(ln,:),'ORIGIN')
    
    data(record_count).Reference{ref_count}.Number = '';
    data(record_count).Reference{ref_count}.Authors = '';
    data(record_count).Reference{ref_count}.Consrtm = '';
    data(record_count).Reference{ref_count}.Title = '';
    data(record_count).Reference{ref_count}.Journal = '';
    data(record_count).Reference{ref_count}.MedLine = '';
    data(record_count).Reference{ref_count}.PubMed = '';
    data(record_count).Reference{ref_count}.Remark = '';
    

    %REFERENCE - Mandatory
    [s,f,t] = regexp(gptext(ln,:),'REFERENCE\s+(\w|\W)+');  %#ok
    if ~isempty(s)
        data(record_count).Reference{ref_count}.Number = strtrim(gptext(ln,t{1}(1):t{1}(2)));
        ln=ln+1;
    end

    %AUTHORS - Mandatory -- though we found examples where this was missing
    [s,f,t] = regexp(gptext(ln,:),'AUTHORS\s+(\w|\W)+');  %#ok
    if ~isempty(s)
        data(record_count).Reference{ref_count}.Authors = strtrim(gptext(ln,t{1}(1):t{1}(2)));
        ln=ln+1;
    end
    
    while ~matchstart(gptext(ln,:),'TITLE') && ~matchstart(gptext(ln,:),'JOURNAL') && ~matchstart(gptext(ln,:),'CONSRTM')
        data(record_count).Reference{ref_count}.Authors=strvcat(data(record_count).Reference{ref_count}.Authors, strtrim(gptext(ln,:))); %#ok
        ln=ln+1;
    end
    
    %CONSRTM - Optional    
    [s,f,t] = regexp(gptext(ln,:),'CONSRTM\s+(\w|\W)+');  %#ok
    if ~isempty(s)
        data(record_count).Reference{ref_count}.Consrtm = strtrim(gptext(ln,t{1}(1):t{1}(2)));
        ln=ln+1;
    end
    while ~matchstart(gptext(ln,:),'TITLE') && ~matchstart(gptext(ln,:),'JOURNAL')
        data(record_count).Reference{ref_count}.Consrtm=strvcat(data(record_count).Reference{ref_count}.Consrtm, strtrim(gptext(ln,:))); %#ok
        ln=ln+1;
    end

    %TITLE - Optional    
    [s,f,t] = regexp(gptext(ln,:),'TITLE\s+(\w|\W)+');  %#ok
    if ~isempty(s)
        data(record_count).Reference{ref_count}.Title = strtrim(gptext(ln,t{1}(1):t{1}(2)));
        ln=ln+1;
    end
    while ~matchstart(gptext(ln,:),'JOURNAL')
        data(record_count).Reference{ref_count}.Title=strvcat(data(record_count).Reference{ref_count}.Title, strtrim(gptext(ln,:))); %#ok
        ln=ln+1;
    end

    %JOURNAL - Mandatory
    [s,f,t] = regexp(gptext(ln,:),'JOURNAL\s+(\w|\W)+');  %#ok
    if ~isempty(s)
        data(record_count).Reference{ref_count}.Journal = strtrim(gptext(ln,t{1}(1):t{1}(2)));
        ln=ln+1;
    end
    
    if matchstart(gptext(ln,:),'REFERENCE')
        ref_count=ref_count+1;
        % next reference
        continue
    end

    if matchstart(gptext(ln,:),'COMMENT') || matchstart(gptext(ln,:),'FEATURES') || matchstart(gptext(ln,:),'BASE COUNT')
        % done with references
        break
    end

    while ~matchstart(gptext(ln,:),'MEDLINE') && ~matchstart(gptext(ln,:),'PUBMED') && ~matchstart(gptext(ln,:),'REMARK')
        data(record_count).Reference{ref_count}.Journal=strvcat(data(record_count).Reference{ref_count}.Journal, strtrim(gptext(ln,:))); %#ok
        ln=ln+1;
        if matchstart(gptext(ln,:),'REFERENCE') || matchstart(gptext(ln,:),'COMMENT') || matchstart(gptext(ln,:),'FEATURES') || matchstart(gptext(ln,:),'BASE COUNT')
            % done with references
            break
        end
    end

    %MEDLINE - Optional    
    [s,f,t] = regexp(gptext(ln,:),'MEDLINE\s+(\d)+');  %#ok
    if ~isempty(s)
        data(record_count).Reference{ref_count}.MedLine = gptext(ln,t{1}(1):t{1}(2));
        ln=ln+1;        
    end

    %PUBMED - Optional    
    [s,f,t] = regexp(gptext(ln,:),'PUBMED\s+(\d+)');  %#ok
    if ~isempty(s)
        data(record_count).Reference{ref_count}.PubMed = gptext(ln,t{1}(1):t{1}(2));
        ln=ln+1;
    end

    %REMARK - Optional    
    [s,f,t] = regexp(gptext(ln,:),'REMARK\s+(\w|\W)+');  %#ok
    if ~isempty(s)
        data(record_count).Reference{ref_count}.Remark = strtrim(gptext(ln,t{1}(1):t{1}(2)));
        ln=ln+1;
    end
    while ~isempty(s) && ~matchstart(gptext(ln,:),'COMMENT') && ~matchstart(gptext(ln,:),'FEATURES') && ~matchstart(gptext(ln,:),'BASE COUNT') && ~matchstart(gptext(ln,:),'REFERENCE')
        data(record_count).Reference{ref_count}.Journal=strvcat(data(record_count).Reference{ref_count}.Journal, strtrim(gptext(ln,t{1}(1):t{1}(2)))); %#ok
        ln=ln+1;
        if matchstart(gptext(ln,:),'COMMENT') || matchstart(gptext(ln,:),'FEATURES') || matchstart(gptext(ln,:),'BASE COUNT'), break, end
    end

    % next reference
    ref_count=ref_count+1;
end

function [ind,flag] = featurelocation(str)
%FEATURELOCATION Translates strings with feature locations into indices.
%
%  [IND,UFB] = FEATURELOCATION(STR) Translates strings with feature
%  locations into indices IND accordingly with the DDBJ/EMBL/GenBank
%  Feature Table Definition (http://www.ncbi.nlm.nih.gov/collab/FT/).
%  Pointers to other records, sites between symbols and inconclusive sites
%  are all parsed as NaNs. Unknown Feature Boundaries are parsed, but UFB
%  will be TRUE if it occurs in at least one section of the descriptor.

% Copyright 2003-2005 The MathWorks, Inc.


if ~isempty(regexp(str,'^complement(','once'))
    [ind,flag] = featurelocation(str(12:end-1));
    ind = fliplr(ind);
elseif ~isempty(regexp(str,'^join(','once'))
    str = str(6:end-1);
    % find subsections of the string (comma separated)
    h = find([1 ~cumsum((str=='(') - (str==')')) & str==',' 1]);
    ind = [];
    flag = false;
    for i = 1:numel(h)-1
        [indT,flagT] = featurelocation(str(h(i):h(i+1)-2));
        ind = [ind indT];
        flag = flag || flagT;
    end
elseif ~isempty(regexp(str,'^.*(?=:)','once')) % points to another record
    ind = [NaN NaN];
    flag = false;
else
    flag = false;
    % sites between symbols cannot be represented with indices
    str = regexprep(str,'\d+\^\d+','NaN');
    % inconclusive locations neither can be represented with indices
    str = regexprep(str,'\(\d+\.\d+\)','NaN');
    % boundaries of unknown feature limits will be parsed, but flagged
    h = regexp(str,'[<>]');
    if ~isempty(h)
        str(h) = '';
        flag = true;
    end
    % is it a range or a single position ?
    if isempty(regexp(str,'\.\.','once'))
        ind = [str2double(str) str2double(str)];
    else
        ind = str2num(strrep(str,'.',' ')); %#ok is a vec
    end
end
    
% REFERENCE: http://www.ncbi.nlm.nih.gov/collab/FT/#3.2
%
% The following is a list of common location descriptors with their meanings: 
% Location                  Description   
% 
% 467                       Points to a single base in the presented sequence 
% 
% 340..565                  Points to a continuous range of bases bounded by and 
%                           including the starting and ending bases
% 
% <345..500                 Indicates that the exact lower boundary point of a 
%                           feature is unknown.  The location begins at some 
%                           base previous to the first base specified (which need 
%                           not be contained in the presented sequence) and con-
%                           tinues to and includes the ending base 
% 
% <1..888                   The feature starts before the first sequenced base 
%                           and continues to and includes base 888
% 
% (102.110)                 Indicates that the exact location is unknown but that 
%                           it is one of the bases between bases 102 and 110, in-
%                           clusive
% 
% (23.45)..600              Specifies that the starting point is one of the bases 
%                           between bases 23 and 45, inclusive, and the end point 
%                           is base 600 
% 
% (122.133)..(204.221)      The feature starts at a base between 122 and 133, 
%                           inclusive, and ends at a base between 204 and 221, 
%                           inclusive
% 
% 123^124                   Points to a site between bases 123 and 124
% 
% 145^177                   Points to a site between two adjacent bases anywhere 
%                           between bases 145 and 177 
% 
% join(12..78,134..202)     Regions 12 to 78 and 134 to 202 should be joined to 
%                           form one contiguous sequence
% 
% complement(join(2691..4571,4918..5163)
%                           Joins regions 2691 to 4571 and 4918 to 5163, then 
%                           complements the joined segments (the feature is 
%                           on the strand complementary to the presented strand)
%  
% join(complement(4918..5163),complement(2691..4571))
%                           Complements regions 4918 to 5163 and 2691 to 4571, 
%                           then joins the complemented segments (the feature is 
%                           on the strand complementary to the presented strand)
%   
% complement(34..(122.126)) Start at one of the bases complementary to those  
%                           between 122 and 126 on the presented strand and finish
%                           at the base complementary to base 34 (the feature is 
%                           on the strand complementary to the presented strand)
% 
% J00194:100..202           Points to bases 100 to 202, inclusive, in the entry 
%                           (in this database) with primary accession number 
%                           'J00194'
%  
