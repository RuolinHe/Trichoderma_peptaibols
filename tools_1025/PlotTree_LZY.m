function PlotTree_LZY(linkm,leaforder,leafnames,inputvisualpara);%,treeorientation,leaftartposi,linewidth,nodenames,treecolor,nodecolors,treelenscale)


hold on;


treevisuapara=[];
treevisuapara.treetowards='Right'; % 'Left','Up','Down'
treevisuapara.start_x_y=[0,0];
treevisuapara.distscale=1;
treevisuapara.treecolor=[0,0,0];
treevisuapara.linewidth=1;
treevisuapara.nodesnamecolor=[0,0,0];
treevisuapara.textshiftsize=1;
treevisuapara.fontsize=15;
if nargin>3 % there is specific visulization argument
    
    
    if isfield(inputvisualpara,'treetowards')
        treevisuapara.treetowards=inputvisualpara.treetowards;
    end
    
    if isfield(inputvisualpara,'start_x_y')
        treevisuapara.start_x_y=inputvisualpara.start_x_y;
    end
    
    if isfield(inputvisualpara,'linewidth')
        treevisuapara.linewidth=inputvisualpara.linewidth;
    end
    if isfield(inputvisualpara,'treecolor')
        treevisuapara.treecolor=inputvisualpara.treecolor;
    end
    if isfield(inputvisualpara,'distscale')
        treevisuapara.distscale=inputvisualpara.distscale;
    end
    if isfield(inputvisualpara,'nodesnamecolor')
        treevisuapara.nodesnamecolor=inputvisualpara.nodesnamecolor;
    end
    
    if isfield(inputvisualpara,'textshiftsize')
        treevisuapara.textshiftsize=inputvisualpara.textshiftsize;
    end
    
else
end



leafnum=size(linkm,1)+1;
nodes_height= zeros(leafnum,1);% along the direction of tree, imaging it as a horizontal tree
nodes_leafposi= zeros(leafnum,1);% along the direction of leafs, imaging it as a horizontal tree

linkage_heights=zeros(leafnum-1,2);% a linkage is described by one heights, two leaf positions
linkage_leafposis=zeros(leafnum-1,2);

for l=1:leafnum
    nodes_leafposi(leaforder(l))=l;
end


for m=1:size(linkm,1) % plot each mother node
    n=m+leafnum;% mother node
    dnodes=linkm(m,1:2);
    
    nodes_height(n)=linkm(m,3);
    nodes_leafposi(n)=mean(nodes_leafposi(dnodes));
    
    % for describing the linkage line
    linkage_heights(m,:)=[nodes_height(n),nodes_height(n)];
    linkage_leafposis(m,:)=nodes_leafposi(dnodes);
end


daughter_heights=zeros(2*leafnum-2,2);
daughter_leafposis=zeros(2*leafnum-2,2);
for d=1:size(daughter_heights,1)
    % location for the mother node
    mloc=find(linkm(:,1)==d|linkm(:,2)==d);
    mid=leafnum+mloc;
    if length(mloc)==1
        daughter_heights(d,:)=[nodes_height(d),nodes_height(mid)];
        daughter_leafposis(d,:)=[nodes_leafposi(d),nodes_leafposi(d)];
    end
end


% plot the tree out
if strcmp(treevisuapara.treetowards,'Left')
    allxs=treevisuapara.distscale*[daughter_heights;linkage_heights];
    allys=[daughter_leafposis;linkage_leafposis];
    
    %set(gca,'Ytick',(allys(1:leafnum)));
    textshift=-1*[treevisuapara.textshiftsize,0];
    
elseif strcmp(treevisuapara.treetowards,'Right')
    allxs=-treevisuapara.distscale*[daughter_heights;linkage_heights];
    allys=[daughter_leafposis;linkage_leafposis];
    
    %set(gca,'Ytick',sort(allys(1:leafnum)),'YAxisLocation','Right');
    textshift=[treevisuapara.textshiftsize,0];
elseif strcmp(treevisuapara.treetowards,'Down')
    allys=treevisuapara.distscale*[daughter_heights;linkage_heights];
    allxs=[daughter_leafposis;linkage_leafposis];
    
    %set(gca,'Ytick',(allys(1:leafnum)));
    textshift=-1*[0,treevisuapara.textshiftsize];
elseif strcmp(treevisuapara.treetowards,'Up')
    allys=-treevisuapara.distscale*[daughter_heights;linkage_heights];
    allxs=[daughter_leafposis;linkage_leafposis];
    
    %set(gca,'Ytick',sort(allys(1:leafnum)),'YAxisLocation','Right');
    textshift=[0,treevisuapara.textshiftsize];
end


allxs=treevisuapara.start_x_y(1)+allxs;
allys=treevisuapara.start_x_y(2)+allys;



for i=1:size(allxs)
    plot(allxs(i,:),allys(i,:),'linewidth',treevisuapara.linewidth,'color',treevisuapara.treecolor)
end

% label the leaf nodes
if ~isempty(leafnames)
    for i=1:leafnum
        th=text(allxs(i)+textshift(1),allys(i)+textshift(2),leafnames{i},'FontSize',treevisuapara.fontsize);
        
        if size(treevisuapara.nodesnamecolor,1)==1
            th.Color=treevisuapara.nodesnamecolor;
        else
            th.Color=treevisuapara.nodesnamecolor(i,:);
        end
        th.HorizontalAlignment ='left ';
        th.FontWeight='bold';
        
    end
end



