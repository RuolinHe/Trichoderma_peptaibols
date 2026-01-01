function hh=Easypcolor_flip_LZY(data,ud,dimgroup,xgroupnames,ygroupnames,dimgrouppro)
% dimgroup should monotonicially increase
% check if it is monotonically increase

%%%%%%%%%%%%%%%%%%%%%
% ud='u': starting from up, ='d': starting from bottom
% dimgroup={[x group id],[y group id]}
% xgroupnames: how that id maps to names
% dimgrouppro: some properties in drawing


hold on;

if ud=='u'
    data=flipud(data);
end
[d_1,d_2]=size(data);
data(d_1+1,:)=nan;
data(:,d_2+1)=nan;

hh=pcolor(data);
set(hh,'LineStyle','none');

% draw bars between groups,
if nargin<3
    dimgroup={[],[]};
end

if nargin<4
    xgroupnames='';
end

if nargin<5
    ygroupnames='';
end



if nargin<6
    dimgrouppro.lw=1;
    dimgrouppro.color=[1,0,0];
end

for dim=1:2
    
    if ud=='u' && dim==1
        lastbarposi=d_1+1;
    else
        lastbarposi=1;
    end
    

    group=dimgroup{dim};
    group_centers=[];
    group_gids=[];
    
    te=[];
    % where the difference start
    for i=2:(1+length(group))
        if i==(1+length(group)) || group(i-1)~=group(i) 
            if dim==1 % first dimension difference, shown in y axis
                
                % position on the x axis, don't change
                d_1nex=1+[0,d_2];
                    
                % position on the y-axis, depends on the direction of y
                if ud=='u' 
                    d_1ney=d_1-[i,i]+2; % from up to down
                    
                else % regular
                    d_1ney=[i,i]; % from down to up
                end
                
                group_centers(end+1)=1/2*(lastbarposi+d_1ney(1));
                
                    
                group_gids(end+1)=group(i-1);
                lastbarposi=d_1ney(1);
                te=[te,lastbarposi];
            else
                d_1nex=[i,i];
                d_1ney=1+[0,d_2];
                
                group_centers(end+1)=1/2*(lastbarposi+d_1nex(1));
                group_gids(end+1)=group(i-1);
                lastbarposi=d_1nex(1);
                te=[te,lastbarposi];
            end
            
            
            plot(d_1nex,d_1ney,'LineWidth',dimgrouppro.lw,'color',dimgrouppro.color)
        end
        
    end % for i=2:length(group)

    [group_centers,group_orders]=sort(group_centers);
    if dim==1 && ~isempty(group_gids) && ~(length(xgroupnames)<max(group_gids))
        set(gca,'YTick',group_centers)
        set(gca,'Yticklabel',xgroupnames(group_gids(group_orders)));
    elseif dim==2 && ~isempty(group_gids) && ~(length(ygroupnames)<max(group_gids))
        set(gca,'XTick',group_centers)
        set(gca,'Xticklabel',ygroupnames(group_gids(group_orders)));
    end
    
end % for dim=1:2

axis([0, d_2, 0, d_1]+1);






