clear
clc

datapath = 'M:\Fusion_dataset_China_Y0120_TIF_Albers\';
outpath = 'N:\2_WheatMap\4_Output\WMI\2001-2020\';

province = {'Shandong','Henan','Hebei','Hubei','Anhui','Jiangsu','Shaanxi','Shanxi','Xinjiang','Gansu','Sichuan'};
[statistic,provincelist] = xlsread('N:\2_WheatMap\2_Input\ProvinceNameArea_2001-2020.xlsx');
years = (2001:2020);

list_prov = dir(datapath);
WWMI = zeros(length(province),length(years));
for ip = 1:length(province)    
    for iq = 1:length(list_prov)
        if strcmp(list_prov(iq).name,province{ip})
            disp(list_prov(iq).name)
            inpath = strcat(datapath,list_prov(iq).name,'\');
            list = dir([inpath,'*.tif']);
            [ndvi_day1,R_ndvi] = geotiffread([inpath,list(1).name]);
            info_ndvi = geotiffinfo([inpath,list(1).name]);
            [row,col,~] = size(ndvi_day1);
            clear ndvi_day1
            
            B1V1_allyear = zeros(length(years),12);
            for iyear = 1:length(years)
                all_ndvi = zeros(row,col,32,'single');
                for i = (46*iyear-9):(46*iyear+22)   % select optimal growing season
                    filename = list(i).name
                    ndvi = geotiffread([inpath,list(i).name]);
                    ndvi = single(ndvi)/1000;
                    ndvi(ndvi>1) = 1;
                    all_ndvi(:,:,i-iyear*46+10) = ndvi;
                end
                
                select1 = zeros(row,col);
                MyPar = parpool('local',40);
                tic
                parfor m = 1:row
                    for n = 1:col
                        grow_series = squeeze(all_ndvi(m,n,:));
                        if sum(grow_series)==0    
                            select1(m,n) = 0;
                        else
                            if  max(grow_series)>0.4
                                select1(m,n) = 1;
                            else
                                select1(m,n) = 0;
                            end
                        end
                    end
                end
                toc
                delete(MyPar);

                select1_series = all_ndvi.* select1;  
                
                select1_max = max(select1_series,[],3); 
                select1_min = min(select1_series,[],3); 
                clear all_ndvi
                
                select1_min1 = select1_min(select1_min~=0);
                select1_max1 = select1_max(select1_max~=0);
                
              
                V1_value1 = prctile(select1_max1,95);
                V1_value2 = prctile(select1_max1,20:20:80);
                V1_value3 = prctile(select1_max1,5);
                V1_value = [V1_value3 V1_value2 V1_value1];
                B1_value1 = prctile(select1_min1,5);
                B1_value2 = prctile(select1_min1,20:20:80);
                B1_value3 = prctile(select1_min1,95);
                B1_value = [B1_value1 B1_value2 B1_value3];
                
                B1V1_allyear(iyear,1:6) = B1_value;
                B1V1_allyear(iyear,7:12) = V1_value;
                
              
                overwinter = select1_series(:,:,26:31);  
                clear select1_series
                aa = [05 20 40 60 80 95];
                bb = [05 20 40 60 80 95];
                
                for kk = 1:length(B1_value)
                    B_value = B1_value(1,kk);
                    for pp = 1:length(V1_value)
                        select2 = zeros(row,col);
                        classified_wheat = zeros(row,col);
                        V_value = V1_value(1,pp);
                        
                        MyPar = parpool('local',40);
                        tic
                        parfor ii = 1:row
                            for jj = 1:col
                                overwinter_series = squeeze(overwinter(ii,jj,:));
                                if sum(overwinter_series) ~= 0
                                    [m1,n1] = max(overwinter_series);
                                    [m2,n2] = min(overwinter_series);
                                    
                                    if (n1<n2)
                                        if m1<V_value
                                            Fv = 1-((V_value - m1)/(V_value-B_value))^2;
                                        else
                                            Fv = 1
                                        end
                                        
                                        if m2>=B_value
                                            Fw = 1-((m2 - B_value)/(V_value-B_value))^2;
                                        else
                                            Fw = 1
                                        end
                                        
                                        D = m1-m2
                                        Fd = 1/(1+exp((V_value-B_value)/2 - D));
                                        
                                        F = Fv*Fw*Fd;
                                        select2(ii,jj) = F;
                                        
                                    else
                                        select2(ii,jj) = -255;
                                    end
                                else
                                    select2(ii,jj) = -255
                                end
                            end
                        end
                        toc
                        delete(MyPar);
                        sum(sum(select2 ~= -255))
                        
                        statisticArea = statistic(ip,iyear);
                        numwheat_Statistic = fix(statisticArea*10000/900);
                        [order,index] = sort(select2(:),'descend');
                        classified_wheat(index(1 : numwheat_Statistic)) = 1;
                        numwheat_classified = sum(classified_wheat(:));
                        classified_wheat = int8(classified_wheat);
                        max(max(classified_wheat))
                        
                        wheat_wwmi = order(numwheat_Statistic,1);
                        WWMI(ip,iyear) = wheat_wwmi;
                        num1 = sprintf('%02d',aa(kk));
                        num2 = sprintf('%02d',bb(pp));
                        
                        outname = strcat(outpath,province{ip},'\',province{ip},'_WheatMap_',num2str(years(iyear)),'_WWMI_Bline',num1,'_Vline',num2,'_DOY289-169_DOY121-161.tif')
                        geotiffwrite(outname,classified_wheat,R_ndvi,'GeoKeyDirectoryTag',info_ndvi.GeoTIFFTags.GeoKeyDirectoryTag);
                        clear  select2 classified_wheat
                        
                    end
                end
            end
            xlswrite([outpath,'1_Bline_Vline_value_2002-2020.xlsx'],B1V1_allyear,province{ip});
        end
        xlswrite([outpath,'AllProvince_2002-2020_WWMI_DOY289-169_DOY97-161.xlsx'],WWMI);
        
        
    end
end

