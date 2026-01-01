function [ai_bi_index,fh]=Silhouette_LZY(groupid,distm,id_2_name,id_2_color)
% remove 'stable' from 'unique' function
unigroups=unique(groupid);
% unigroups=unique(groupid,'stable');
members_eachgroup=[];
for i=1:length(unigroups)
    members_eachgroup{i}=find(groupid==unigroups(i)); % it needs loc as input
end

for i=1:length(unigroups)
    thisgroupid=unigroups(i);
    thisgroupmembers=members_eachgroup{i};
    
    for j=1:length(thisgroupmembers)
        member=thisgroupmembers(j);
        samegroup_members=setdiff(thisgroupmembers,member);
        
        % ai
         if ~isempty(samegroup_members)
            ai_bi_index(member,1)=mean(distm(member,samegroup_members)','omitnan'); % within the same
         end
        
         % bi
         meandist_withothergroup=[];
          for o=setdiff(1:length(unigroups),i)
              othergroup_members=members_eachgroup{o};
              meandist_withothergroup(end+1)=mean(distm(member,othergroup_members)','omitnan'); % within the same
          end
         ai_bi_index(member,2)=min(meandist_withothergroup);
    end
end

ai=ai_bi_index(:,1);
bi=ai_bi_index(:,2);
ai_bi_index(:,3)=(bi-ai)./max([ai,bi],[],2);


% want to plot it out
if nargin>2
    
    fignum = length(findobj('type','figure'));
    
    fh=figure(fignum+1);
    hold on;
    
    y_low=0;
    ymiddle=zeros(1,length(unigroups));
    for i=1:length(unigroups)
        gid=unigroups(i);
        members=find(groupid==gid);

            barh(y_low-(1:length(members)),sort(ai_bi_index(members,3)),'EdgeColor','none','FaceColor',id_2_color(gid,:));
        ymiddle(1,i)=y_low-1/2*length(members);
        y_low=y_low-length(members)-1;
        
    end
    set(gca,'ytick',fliplr(ymiddle),'yticklabel',id_2_name(unigroups(length(unigroups):-1:1)));
    
    figure;
end




