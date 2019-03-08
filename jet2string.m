clc; clear all;close all;
path_adt = '..\allsat_sla\';
folder_list=[1993:1:2015];
lat_acc = ncread('./ofes_jet_detection_SO_2015.nc','lat'); 
lon_acc = ncread('./ofes_jet_detection_SO_2015.nc','lon');
for y = 1:length(folder_list)
    jet_acc = ncread(['./ofes_jet_detection_SO_',num2str(folder_list(y)),'.nc'],'jet_locations');
    jet_year = [];
    data_year=[];

    for d = 1:size(jet_acc,3)
        d+y*1000
        jet_cache = squeeze(jet_acc(:,:,d));
        data=[];
        [xx,yy] = find(jet_cache==1);
        data = [xx,yy];
        [data(:,1),ind] = sort(data(:,1));
        data(:,2) = data(ind,2);
        it_1 = find(data(:,1)==1);
        it_1440 = find(data(:,1)==1440);
        
        distmat=zeros(size(data,1),size(data,1));
        for i=1:size(data,1)
            for j=i:size(data,1)
                distmat(i,j)=sqrt((data(i,1:2)-data(j,1:2))*(data(i,1:2)-data(j,1:2))');
            end
        end
        for i = min(it_1):max(it_1)
            for j = min(it_1440):max(it_1440)
                distmat(i,j)=sqrt((data(i,1:2)+[1440,0]-data(j,1:2))*(data(i,1:2)+[1440,0]-data(j,1:2))');
            end
        end
        for i=1:size(data,1)
            for j=i:size(data,1)
                distmat(j,i)=distmat(i,j);
            end
        end
        distmat(distmat==0)=nan;
        
        jet_index = zeros(length(data),1);        
        for k = 1:length(data)
            it = find(distmat(k,k+1:end)<=1.6);  %%% 距离阈值可能过小，跨360经度线连不上
            if length(it)>0
                if jet_index(k)==0
                    jet_index(k) = k;
                    jet_index(it+k) = jet_index(k);
                else
                    jet_index(it+k) = jet_index(k);
                end
            end
        end
        
        data_year{d} = data;
        jet_year{d} = jet_index;
    end
    save(['jet_long',num2str(folder_list(y)),'.mat'],'jet_year','data_year');
end

%% 修正跨360度经度线
clc;clear all; close all;
% dir_mat = dir('./*.mat');
for y = 1:length(folder_list)
    load(['jet_long',num2str(y+1992),'.mat']);
    jet_long_year = [];
    for d = 1:length(jet_year)
        num_initial = unique(jet_year{d});
        jet_length=[];
        jet_long=[];
        data = data_year{d};
        for n = 2:length(num_initial)    % num_initial==0 不是成串的急流
            jet_length(n-1) = length(find(jet_year{d}==num_initial(n)));
        end
        h=1;k=1;it_1=[];it_1440=[];
        for i = 2:length(num_initial)
            it = find(jet_year{d} == num_initial(i));
            jet_long{i-1,1} = data(it,:);
            if length(find(data(it,1)==1)) > 0
                it_1(h)=i-1;
                h=h+1;
            end
            if length(find(data(it,1)==1440)) > 0
                it_1440(k)=i-1;
                k=k+1;
            end
        end
        jet_long1 = jet_long(it_1);
        jet_long1440 = jet_long(it_1440);
        dist=[];
        for v = 1:length(it_1440)
            cache1440 = jet_long1440{v};
            for w = 1:length(it_1)
                cache1 = jet_long1{w};
            dist(v,w) = sqrt((cache1440(cache1440(:,1)==1440,:) - cache1(cache1(:,1)==1,:)-[1440,0])...
                *(cache1440(cache1440(:,1)==1440,:) - cache1(cache1(:,1)==1,:)-[1440,0])');
            end
        end
        [a,b] = find(dist<=1.6);
%         it_1440(a);
%         it_1(b);
        jet_long_revise = jet_long;
        for l = 1:length(b)
            jet_long_revise{it_1(b(l))} = [jet_long{it_1440(a(l))}; jet_long{it_1(b(l))}];
        end
        jet_long_revise(it_1440(a))=[];
        
%         for i = 2:length(num_initial)
%             it = find(jet_year{d} == num_initial(i));
%             jet_long{i-1,1} = data(it,:);
%         end
        jet_long_year{d,1} = jet_long_revise;
    end
    save(['jet_long_revise',num2str(y+1992),'.mat'],'jet_long_year')
end

%%

% data = data_year{2};
% abc = unique(jet_year{2});
% jet_length=[];
% for n = 1:length(abc)
%     jet_length(n) = length(find(jet_year{2}==abc(n)));
% end
% 
% figure
% m_proj('stereographic','lat',-90,'long',0,'radius',60);
% m_grid('xtick',6,'tickdir','out','ytick',[-85 -65],'linest','-.');
% m_coast('patch',[255 251 71]/255,'edgecolor','k');
% % m_proj('miller','lat',[-85.5 -30],'long',[0.5 359.5]);
% % m_proj('miller','lat',[-90 -30],'long',[0 60]);
% % m_grid('xtick',9,'ytick',9,'tickdir','out','linest',':');
% hold on
% for n = 2:length(abc)
%     if jet_length(n)>=2
%         lat_test = lat_acc(data(find(jet_index==abc(n)),2));
%         lon_test = lon_acc(data(find(jet_index==abc(n)),1));
%         [lon_test,ind] = sort(lon_test);
%         lat_test = lat_test(ind);
%         m_plot(lon_test,lat_test,'-','linewidth',1)
%     end
% end
% print('-dtiff','-r600','nextday')
% 
% 
%  plot(data(:,2),data(:,1),'bo','markersize',2)  
%  hold on
% for n = 2:length(abc)
%     if jet_length(n)>=2
%         lat_test = data(find(jet_index==abc(n)),2);
%         lon_test = data(find(jet_index==abc(n)),1);
%         [lon_test,ind] = sort(lon_test);
%         lat_test = lat_test(ind);
%         plot(lat_test,lon_test,'-','linewidth',1)
%     end
% end