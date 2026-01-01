function [fhandles,optleafOrder,ai_bi_index]=DistMatrix_Analysis_HRL20220210(distm,figoutput_list,groups,figtitle)
% in "figoutput_list",
% "hist" plot the histogram of distance;
% 'heatmap_general' plot the structural in general
% 'phytree_single' plot the phylogenetic tree by single linkage (default)
% 'phytree_average' plot the phylogentic tree by average linkage
% groups: is different strings about the grouping method
%
% updated by HRL. Now 'groups' could have a filed named 'matrix', and show
% mixed properties such as substrate promiscuity
% remove 'stable' from 'unique' function
fignum = length(findobj('type','figure'));

subplotnum=[50,10,5];

if nargin<2
    figoutput_list={'hist','heatmap_general','phytree_single','Silhouette'};
end
if nargin<3
    groupmethods=0;
else
    groupmethods=length(groups);
end
if nargin<4
    figtitle='';
end

matrixsize=size(distm,1);
if distm(1,1)~=0
    for i=1:matrixsize
        distm(i,i)=0;
    end
end

fhandles=[];
optleafOrder=1:size(distm,1);

distline=squareform(distm);

% plot the histogram
if sum(strcmp('hist',figoutput_list))>0
    histnum=max(10,round(length(distline)/5));
    fignum=fignum+1;fhandles(end+1)=figure(fignum);
    hist(distline,histnum);
    title(figtitle)
    xlabel('Distance');
end

% plot the overall distance heatmap
if sum(strcmp('heatmap_general',figoutput_list))>0
    fignum=fignum+1;fhandles(end+1)=figure(fignum);
    Easypcolor_flip_LZY(distm,'u');
    colorbar;
    title(figtitle)
    caxis([min(distline),max(distline)]);
end

% plot the phytree and corresponding ordering
phytreeid=0;
phytreelocs=find(strncmp('phytree',figoutput_list,length('phytree')));
for i=1:length(phytreelocs)
    phystr=figoutput_list{phytreelocs};
    strlist=strsplit(phystr,'_');
    if length(strlist)==2
        treemethod=strlist{2};
    else
        treemethod='single';
    end
    
    distmforphy=distm;
    distmforphy(isnan(distmforphy))=mean(distm(:),'omitnan');
    
    % get the phylogenetic tree of the distance matrix
    linkm = linkage(distmforphy,treemethod);
    optleafOrder = optimalleaforder(linkm,distmforphy);
    
    % plot the heatmap
    fignum=fignum+1;fhandles(end+1)=figure(fignum);phytreeid=fignum; hold on;
    
    subplot(1,subplotnum(1),1:subplotnum(2))
    
    treevisuapara=[];
    treevisuapara.start_x_y=[0,0.5];
    PlotTree_LZY(linkm,fliplr(optleafOrder),'',treevisuapara)
    ylim([0,length(optleafOrder)+1])
    axis off
    
    subplot(1,subplotnum(1),(subplotnum(2)+groupmethods*subplotnum(3)+1):(subplotnum(1)))
    Easypcolor_flip_LZY(distm(optleafOrder,optleafOrder),'u');
    ch=colorbar;
    ylabel(ch,[figtitle,' Dist'])
    caxis([min(distline),max(distline)]);
    ylim([0,length(optleafOrder)+1])
    axis off
    
end

% how the distance ~ groups;

fignum=fignum+1;fhandles(end+1)=figure(fignum);% this is for the legend
legendid=fignum;

for gm=1:groupmethods
    grouppro=groups{gm};
    group_ids=grouppro.ids; % should be number
    id_2_name=grouppro.id_2_name;
    id_2_color=grouppro.id_2_color;
    group_matrix=[];
    if isfield(grouppro,'matrix')
        group_matrix=grouppro.matrix;
    end
    %  group_identifiers=groups{gm};
    [unigroup_identifiers,~,id_in_uni]=unique(group_ids); % always use id_in_uni as identification for groups
%     [unigroup_identifiers,~,id_in_uni]=unique(group_ids,'stable'); % always use id_in_uni as identification for groups
    
    %%%%%%%%%%% output color legends  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(id_2_color)
        figure(legendid)
        subplot(1,groupmethods,gm)
        hold on;
        for i=1:length(id_2_color)
            if isempty(group_matrix)
                id=unigroup_identifiers(i);
            else
                id = i;
            end
            ph=patch([0,1,1,0],[0,0,1,1]*0.8-i,id_2_color(id,:));
            text(1.1,0.5-i,id_2_name{id})
            ph.EdgeColor='none';
        end
        axis equal
        title(grouppro.name)
        axis off
    end
    
    %%%%%%%%%%% color the phytree if it exist  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(optleafOrder) && ~isempty(id_2_color)
        figure(phytreeid)
        subplot(1,subplotnum(1),(subplotnum(2)+1+(gm-1)*subplotnum(3)):(subplotnum(2)+gm*subplotnum(3)-1))
        hold on;
        if ~isempty(group_matrix) % if it is a matrix
            for i=1:length(optleafOrder)
                id=optleafOrder(i);
                loc_index=find(group_matrix(id,:));
                loc_sum=0;
                for ii=1:length(loc_index)
                    new_sum=loc_sum+group_matrix(id,loc_index(ii));
                    ph=patch([loc_sum,new_sum,new_sum,loc_sum],length(optleafOrder)+1+[0,0,1,1]-i,id_2_color(loc_index(ii),:));
                    ph.EdgeColor='none';
                    loc_sum=new_sum;
                end
            end
        else
            for i=1:length(optleafOrder)
                id=group_ids(optleafOrder(i));
                ph=patch([0,1,1,0],length(optleafOrder)+1+[0,0,1,1]-i,id_2_color(id,:));
                ph.EdgeColor='none';
            end
        end
        axis off
        xlim([0,1])
        ylim([0,length(optleafOrder)+1])
        title(grouppro.name)
    end
    
    %%%%%%%%%%% re-order the heatmap by groups, then plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    neworder=[];
    neworder_gid=[];
    
    locinleaf=zeros(size(distm));locinleaf(optleafOrder)=1:length(optleafOrder);
    groupsize=zeros(1,length(unigroup_identifiers));
    accumusize=zeros(1,length(unigroup_identifiers));
    for i=1:length(unigroup_identifiers)
        toadd=find(id_in_uni==i);
        
        groupsize(i)=length(toadd);
        if i>1
            accumusize(i)=accumusize(i-1)+groupsize(i-1);
        end
        % re-order these "to add" by the optimal leaf
        %[~,order]=sort(locinleaf(toadd));
        %toadd=toadd(order);
        neworder=[neworder;toadd];
        neworder_gid=[neworder_gid;unigroup_identifiers(i)*ones(groupsize(i),1)];
        
    end
    
    fignum=fignum+1;
    fhandles(end+1)=figure(fignum);
    wholewindow=20;
    if ~isempty(id_2_color)
        largewindow=19;
    else
        largewindow=wholewindow;
    end
    smallwindow=(wholewindow-largewindow);
    allsubplots=reshape(1:(wholewindow^2),wholewindow,wholewindow)';
    upplots=allsubplots(1:smallwindow,1:largewindow);
    rightwindow=allsubplots((smallwindow+1):wholewindow,(largewindow+1):wholewindow);
    mainplots=allsubplots((smallwindow+1):wholewindow,1:largewindow);
    
    subplot(wholewindow,wholewindow, mainplots(:)'  );
    Easypcolor_flip_LZY(distm(neworder,neworder),'u',{neworder_gid,neworder_gid},id_2_name,id_2_name);
    xlim([1,size(distm,1)+1])
    caxis([min(distline),max(distline)]);

    ch=colorbar;
    ylabel(ch,[figtitle,' Dist'])
    
    title(grouppro.name)

    if ~isempty(id_2_color)
        alllen=size(distm,1)+1;
        subplot(wholewindow,wholewindow,rightwindow(:)')
        hold on;
         for i=1:length(unigroup_identifiers)
             id=unigroup_identifiers(i);
             ph=patch([0,0,1,1],[1,1,1,1]*alllen-[0,1,1,0]*groupsize(i)-accumusize(i),id_2_color(id,:));
            ph.EdgeColor='none';
        end
        xlim([0 1])
        ylim([1,alllen])
        axis off

        subplot(wholewindow,wholewindow, upplots(:)')
        hold on;
        for i=1:length(unigroup_identifiers)
            id=unigroup_identifiers(i);
            ph=patch(1+[0,1,1,0]*groupsize(i)+accumusize(i),[0,0,1,1],id_2_color(id,:));
            ph.EdgeColor='none';
        end
        ylim([0 1])
        xlim([1,alllen])
        axis off
    end
    
    %%%% the average similarity between groups
     meandist_bygroup=nan*zeros(length(unigroup_identifiers));
     for g1=1:length(unigroup_identifiers)
         members_g1=find(id_in_uni==g1);
         for g2=1:g1
             members_g2=find(id_in_uni==g2);
             if g1==g2 && length(members_g1)>1 && length(members_g2)>1
                 dist_in_line=squareform(distm(members_g1,members_g2));
             else
                 dist_in_line=distm(members_g1,members_g2);dist_in_line=dist_in_line(:);
             end
             meandist_bygroup(g1,g2)=mean(dist_in_line,'omitnan');
             meandist_bygroup(g2,g1)= meandist_bygroup(g1,g2);
         end
     end
      
    fignum=fignum+1;
    fhandles(end+1)=figure(fignum);
    
    dimgrouppro=[];
    dimgrouppro.lw=1;
    dimgrouppro.color=[1,1,1];
    
    Easypcolor_flip_LZY(meandist_bygroup,'u',{unigroup_identifiers,unigroup_identifiers},id_2_name,id_2_name,dimgrouppro);
    title([grouppro.name,' meandist'])
    caxis([min(distline),max(distline)]);
    ch=colorbar;
    
     %%%%%%%%%%% re-order the heatmap by groups, then plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %need to use neworder, and neworder_gid=[];
     if sum(strcmp('Silhouette',figoutput_list))>0
         [ai_bi_index,fh]=Silhouette_LZY(group_ids,distm,id_2_name,id_2_color);
         fignum=fignum+1;
    end
   
    
end % for gm=1:groupmethods


