sepstr='/';
toolboxpath='../tools_1025/';
if exist(toolboxpath,'dir')
    addpath(toolboxpath)
else
    error('Wrong path for MATLAB additional toolbox, please check\n')
end
if ~exist('best_peptaibol','var')
    load("Known_peptaibol_struct.mat")
end
%% A4-A5
Tex2_index=find(ismember(best_peptaibol.peptaibol_name,'TEX2'));
groups=[];
for i=1:4
    grouppro=[];
     if i==1
        grouppro.name='Specificity';
        substrate_mode=best_peptaibol.substrate_mode;
        substrate_mode(best_peptaibol.peptaibol_list==Tex2_index)=[];% remove short peptaibol
        specificity_freq_matrix=best_peptaibol.specificity_freq_matrix;
        specificity_freq_matrix(best_peptaibol.peptaibol_list==Tex2_index,:)=[];
        loc_specificity_freq_matrix_sum=sum(specificity_freq_matrix);
        unique_substrate=best_peptaibol.unique_substrate;
        empty_substrate=find(loc_specificity_freq_matrix_sum==0);
        unique_substrate(empty_substrate)=[];
        specificity_freq_matrix(:,empty_substrate)=[];
        add_value=zeros(size(substrate_mode));
        for j = 1:length(empty_substrate)
            add_value(substrate_mode>empty_substrate(j))=add_value(substrate_mode>empty_substrate(j))-1;
        end
        substrate_mode=substrate_mode+add_value;
        [uni_substrate_mode,~,grouppro.ids]=unique(substrate_mode);
        non_mode_substrate=setdiff(unique(grouppro.ids),uni_substrate_mode);
        specificity_freq_matrix=[specificity_freq_matrix,specificity_freq_matrix(:,non_mode_substrate)];
        specificity_freq_matrix(:,non_mode_substrate)=[];
        unique_substrate=[unique_substrate;unique_substrate(non_mode_substrate)];
        unique_substrate(non_mode_substrate)=[];
        grouppro.matrix=specificity_freq_matrix;
        grouppro.id_2_name=unique_substrate;
        hex_str={'#FD5D67','#FFB366','#7EE1A8','#CFE9B0','#B8E001','#3C7747','#AA7A38','#6C5BFF','#D3C1B3','#00CCCC','#333333','#7A5D8B','#818689'};
        hex_str=[hex_str,hex_str(non_mode_substrate)];
        hex_str(non_mode_substrate)=[];
        grouppro.id_2_color=hex2rgb(hex_str);
    elseif i==2
        module_ancestor_list=best_peptaibol.module_ancestor_list;
        module_ancestor_list(best_peptaibol.peptaibol_list==Tex2_index)=[];
        [grouppro.id_2_name,~,grouppro.ids] = unique(module_ancestor_list);
        grouppro.name='Ancestor Module';
        grouppro.id_2_name=cellstr(num2str(grouppro.id_2_name));
        hex_str={'#FFD9B3','#FFBF80','#FF8000','#D7261C','#990000','#4B4FC5','#6A63D7','#6BC1E9','#5AAAE7','#4A90E2','#3D7AD6','#3264C8','#2A55B5','#D9F2B3','#BFE680','#8CD790','#B8E001','#00B300','#008080','#3C7747','#FFF2B3','#CCAA80'};
        grouppro.id_2_color=hex2rgb(hex_str);
    elseif i==3
        module_list=best_peptaibol.module_list;
        module_list(best_peptaibol.peptaibol_list==Tex2_index)=[];
        [grouppro.id_2_name,~,grouppro.ids] = unique(module_list);
        grouppro.name='Module';
        grouppro.id_2_name=cellstr(num2str(grouppro.id_2_name));
        hex_str={'#FFD9B3','#FFBF80','#FF8000','#D7261C','#990000','#4B4FC5','#6A63D7','#6BC1E9','#5AAAE7','#4A90E2','#3D7AD6','#3264C8','#2A55B5','#D9F2B3','#BFE680','#8CD790','#B8E001','#00B300','#008080','#3C7747'};
        grouppro.id_2_color=hex2rgb(hex_str);
    elseif i==4
        peptaibol_list=best_peptaibol.peptaibol_list;
        peptaibol_list(best_peptaibol.peptaibol_list==Tex2_index)=[];
        peptaibol_list(peptaibol_list>Tex2_index)=peptaibol_list(peptaibol_list>Tex2_index)-1;
        peptaibol_name=best_peptaibol.peptaibol_name_str;
        peptaibol_name(Tex2_index)=[];
        grouppro.name='PS';
        grouppro.ids=peptaibol_list;
        grouppro.id_2_name=peptaibol_name;
        hex_str={'#FFBF80','#D9F2B3','#00B300','#B3D9FF','#BFE680','#FF8000','#883890','#008FB3','#FFD9B3','#0066CC','#990000','#CCAA80','#008080','#FFF2B3','#9980FF','#6600CC','#D9CCFF','#D7261C'};
        grouppro.id_2_color=hex2rgb(hex_str);
    end
    groups{i}=grouppro;
end
A4_A5_dist_matrix=best_peptaibol.A4_A5_dist_matrix;
A4_A5_dist_matrix(best_peptaibol.peptaibol_list==Tex2_index,:)=[];
A4_A5_dist_matrix(:,best_peptaibol.peptaibol_list==Tex2_index)=[]; 
[~,optleafOrder,~]=DistMatrix_Analysis_HRL20220210(A4_A5_dist_matrix,{'phytree','Silhouette'},groups,'A4-A5 (p-distance)');
%% 34AA
Tex2_index=find(ismember(best_peptaibol.peptaibol_name,'TEX2'));
groups=[];
for i=1:4
    grouppro=[];
     if i==1
        grouppro.name='Specificity';
        substrate_mode=best_peptaibol.substrate_mode;
        substrate_mode(best_peptaibol.peptaibol_list==Tex2_index)=[];% remove short peptaibol
        specificity_freq_matrix=best_peptaibol.specificity_freq_matrix;
        specificity_freq_matrix(best_peptaibol.peptaibol_list==Tex2_index,:)=[];
        loc_specificity_freq_matrix_sum=sum(specificity_freq_matrix);
        unique_substrate=best_peptaibol.unique_substrate;
        empty_substrate=find(loc_specificity_freq_matrix_sum==0);
        unique_substrate(empty_substrate)=[];
        specificity_freq_matrix(:,empty_substrate)=[];
        add_value=zeros(size(substrate_mode));
        for j = 1:length(empty_substrate)
            add_value(substrate_mode>empty_substrate(j))=add_value(substrate_mode>empty_substrate(j))-1;
        end
        substrate_mode=substrate_mode+add_value;
        [uni_substrate_mode,~,grouppro.ids]=unique(substrate_mode);
        non_mode_substrate=setdiff(unique(grouppro.ids),uni_substrate_mode);
        specificity_freq_matrix=[specificity_freq_matrix,specificity_freq_matrix(:,non_mode_substrate)];
        specificity_freq_matrix(:,non_mode_substrate)=[];
        unique_substrate=[unique_substrate;unique_substrate(non_mode_substrate)];
        unique_substrate(non_mode_substrate)=[];
        grouppro.matrix=specificity_freq_matrix;
        grouppro.id_2_name=unique_substrate;
        hex_str={'#FD5D67','#FFB366','#7EE1A8','#CFE9B0','#B8E001','#3C7747','#AA7A38','#6C5BFF','#D3C1B3','#00CCCC','#333333','#7A5D8B','#818689'};
        hex_str=[hex_str,hex_str(non_mode_substrate)];
        hex_str(non_mode_substrate)=[];
        grouppro.id_2_color=hex2rgb(hex_str);
    elseif i==2
        module_ancestor_list=best_peptaibol.module_ancestor_list;
        module_ancestor_list(best_peptaibol.peptaibol_list==Tex2_index)=[];
        [grouppro.id_2_name,~,grouppro.ids] = unique(module_ancestor_list);
        grouppro.name='Ancestor Module';
        grouppro.id_2_name=cellstr(num2str(grouppro.id_2_name));
        hex_str={'#FFD9B3','#FFBF80','#FF8000','#D7261C','#990000','#4B4FC5','#6A63D7','#6BC1E9','#5AAAE7','#4A90E2','#3D7AD6','#3264C8','#2A55B5','#D9F2B3','#BFE680','#8CD790','#B8E001','#00B300','#008080','#3C7747','#FFF2B3','#CCAA80'};
        grouppro.id_2_color=hex2rgb(hex_str);
    elseif i==3
        module_list=best_peptaibol.module_list;
        module_list(best_peptaibol.peptaibol_list==Tex2_index)=[];
        [grouppro.id_2_name,~,grouppro.ids] = unique(module_list);
        grouppro.name='Module';
        grouppro.id_2_name=cellstr(num2str(grouppro.id_2_name));
        hex_str={'#FFD9B3','#FFBF80','#FF8000','#D7261C','#990000','#4B4FC5','#6A63D7','#6BC1E9','#5AAAE7','#4A90E2','#3D7AD6','#3264C8','#2A55B5','#D9F2B3','#BFE680','#8CD790','#B8E001','#00B300','#008080','#3C7747'};
        grouppro.id_2_color=hex2rgb(hex_str);
    elseif i==4
        peptaibol_list=best_peptaibol.peptaibol_list;
        peptaibol_list(best_peptaibol.peptaibol_list==Tex2_index)=[];
        peptaibol_list(peptaibol_list>Tex2_index)=peptaibol_list(peptaibol_list>Tex2_index)-1;
        peptaibol_name=best_peptaibol.peptaibol_name_str;
        peptaibol_name(Tex2_index)=[];
        grouppro.name='PS';
        grouppro.ids=peptaibol_list;
        grouppro.id_2_name=peptaibol_name;
        hex_str={'#FFBF80','#D9F2B3','#00B300','#B3D9FF','#BFE680','#FF8000','#883890','#008FB3','#FFD9B3','#0066CC','#990000','#CCAA80','#008080','#FFF2B3','#9980FF','#6600CC','#D9CCFF','#D7261C'};
        grouppro.id_2_color=hex2rgb(hex_str);
    end
    groups{i}=grouppro;
end
A4_A5_dist_matrix=best_peptaibol.AA34_p_distance_dist;
A4_A5_dist_matrix(best_peptaibol.peptaibol_list==Tex2_index,:)=[];
A4_A5_dist_matrix(:,best_peptaibol.peptaibol_list==Tex2_index)=[]; 
[~,optleafOrder,~]=DistMatrix_Analysis_HRL20220210(A4_A5_dist_matrix,{'phytree','Silhouette'},groups,'34AA (p-distance)');
%% unknown version
Tex2_index=find(ismember(best_peptaibol.peptaibol_name,'TEX2'));
hypoxylon_index=find(ismember(best_peptaibol.Strain_used,'Trichoderma hypoxylon'));
groups=[];
for i=1:4
    grouppro=[];
     if i==1
        grouppro.name='Specificity';
        substrate_mode=best_peptaibol.substrate_mode;
        substrate_mode(best_peptaibol.peptaibol_list==hypoxylon_index)=length(best_peptaibol.unique_substrate)+1;%Unknown
        substrate_mode(best_peptaibol.peptaibol_list==Tex2_index)=[];% remove short peptaibol
        specificity_freq_matrix=best_peptaibol.specificity_freq_matrix;
        specificity_freq_matrix(best_peptaibol.peptaibol_list==hypoxylon_index,:)=0;
        specificity_freq_matrix(best_peptaibol.peptaibol_list==Tex2_index,:)=[];
        loc_specificity_freq_matrix_sum=sum(specificity_freq_matrix);
        unique_substrate=best_peptaibol.unique_substrate;
        empty_substrate=find(loc_specificity_freq_matrix_sum==0);
        unique_substrate(empty_substrate)=[];
        specificity_freq_matrix(:,empty_substrate)=[];
        specificity_freq_matrix=[specificity_freq_matrix,[zeros(length(best_peptaibol.substrate_mode)-sum(best_peptaibol.peptaibol_list==hypoxylon_index)-sum(best_peptaibol.peptaibol_list==Tex2_index),1);ones(sum(best_peptaibol.peptaibol_list==hypoxylon_index),1)]];
        unique_substrate=[unique_substrate;{'Unknown'}];
        add_value=zeros(size(substrate_mode));
        for j = 1:length(empty_substrate)
            add_value(substrate_mode>empty_substrate(j))=add_value(substrate_mode>empty_substrate(j))-1;
        end
        substrate_mode=substrate_mode+add_value;
        [uni_substrate_mode,~,grouppro.ids]=unique(substrate_mode);
        non_mode_substrate=setdiff(unique(grouppro.ids),uni_substrate_mode);
        specificity_freq_matrix=[specificity_freq_matrix,specificity_freq_matrix(:,non_mode_substrate)];
        specificity_freq_matrix(:,non_mode_substrate)=[];
        unique_substrate=[unique_substrate;unique_substrate(non_mode_substrate)];
        unique_substrate(non_mode_substrate)=[];
        grouppro.matrix=specificity_freq_matrix;
        grouppro.id_2_name=unique_substrate;
        hex_str={'#FD5D67','#FFB366','#7EE1A8','#CFE9B0','#B8E001','#3C7747','#AA7A38','#6C5BFF','#D3C1B3','#00CCCC','#7A5D8B','#818689','#FFFFFF'};% remove Tyr and add Unknown
        hex_str=[hex_str,hex_str(non_mode_substrate)];
        hex_str(non_mode_substrate)=[];
        grouppro.id_2_color=hex2rgb(hex_str);
    elseif i==2
        module_ancestor_list=best_peptaibol.module_ancestor_list;
        module_ancestor_list(best_peptaibol.peptaibol_list==Tex2_index)=[];
        [grouppro.id_2_name,~,grouppro.ids] = unique(module_ancestor_list);
        grouppro.name='Ancestor Module';
        grouppro.id_2_name=cellstr(num2str(grouppro.id_2_name));
        hex_str={'#FFD9B3','#FFBF80','#FF8000','#D7261C','#990000','#4B4FC5','#6A63D7','#6BC1E9','#5AAAE7','#4A90E2','#3D7AD6','#3264C8','#2A55B5','#D9F2B3','#BFE680','#8CD790','#B8E001','#00B300','#008080','#3C7747','#FFF2B3','#CCAA80'};
        grouppro.id_2_color=hex2rgb(hex_str);
    elseif i==3
        module_list=best_peptaibol.module_list;
        module_list(best_peptaibol.peptaibol_list==Tex2_index)=[];
        [grouppro.id_2_name,~,grouppro.ids] = unique(module_list);
        grouppro.name='Module';
        grouppro.id_2_name=cellstr(num2str(grouppro.id_2_name));
        hex_str={'#FFD9B3','#FFBF80','#FF8000','#D7261C','#990000','#4B4FC5','#6A63D7','#6BC1E9','#5AAAE7','#4A90E2','#3D7AD6','#3264C8','#2A55B5','#D9F2B3','#BFE680','#8CD790','#B8E001','#00B300','#008080','#3C7747'};
        grouppro.id_2_color=hex2rgb(hex_str);
    elseif i==4
        peptaibol_list=best_peptaibol.peptaibol_list;
        peptaibol_list(best_peptaibol.peptaibol_list==Tex2_index)=[];
        peptaibol_list(peptaibol_list>Tex2_index)=peptaibol_list(peptaibol_list>Tex2_index)-1;
        peptaibol_name=best_peptaibol.peptaibol_name_str;
        peptaibol_name(Tex2_index)=[];
        grouppro.name='PS';
        grouppro.ids=peptaibol_list;
        grouppro.id_2_name=peptaibol_name;
        hex_str={'#FFBF80','#D9F2B3','#00B300','#B3D9FF','#BFE680','#FF8000','#883890','#008FB3','#FFD9B3','#0066CC','#990000','#CCAA80','#008080','#FFF2B3','#9980FF','#6600CC','#D9CCFF','#D7261C'};
        grouppro.id_2_color=hex2rgb(hex_str);
    end
    groups{i}=grouppro;
end
A4_A5_dist_matrix=best_peptaibol.A4_A5_dist_matrix;
A4_A5_dist_matrix(best_peptaibol.peptaibol_list==Tex2_index,:)=[];
A4_A5_dist_matrix(:,best_peptaibol.peptaibol_list==Tex2_index)=[]; 
[~,optleafOrder,~]=DistMatrix_Analysis_HRL20220210(A4_A5_dist_matrix,{'phytree','Silhouette'},groups,'A4-A5 (p-distance)');
%% predict unknown substrates
unknown_id=find(ismember(groups{1}.id_2_name,'Unknown'));
unknown_index=find(groups{1}.matrix(:,unknown_id));
closed_substrate_matrix=zeros(length(unknown_index),size(groups{1}.matrix,2));
closed_dist_list=zeros(length(unknown_index),1);
for i = 1:length(unknown_index)
    loc_dist=A4_A5_dist_matrix(:,unknown_index(i));
    [~,I]=sort(loc_dist);
    for j = 2:length(I) % 最近的是它自己
        if ~ismember(I(j),unknown_index)
            closed_substrate_matrix(i,:)=groups{1}.matrix(I(j),:);
            closed_dist_list(i)=A4_A5_dist_matrix(i,I(j));
            break
        end
    end
end
%%
predicted_substrate=cell(length(unknown_index),1);
for i = 1:length(unknown_index)
    predicted_substrate(i)=join(groups{1}.id_2_name(closed_substrate_matrix(i,:)>0),'/');
end
%%
tmp_sheet=[{'Module index'},groups{1}.id_2_name',{'Closest distance','Predicted substrate'};num2cell([(1:length(unknown_index))',closed_substrate_matrix,closed_dist_list]),predicted_substrate];
writecell(tmp_sheet,'Prediction_unknown.xlsx')