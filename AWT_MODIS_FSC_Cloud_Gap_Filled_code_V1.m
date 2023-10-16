%% ********************************************************************* %%
% AWT MODIS FSC Cloud Gap Filled product code
% 0-100:FSC, 237:Water£¬250:Cloud£¬253/255:Nodata
% Code by Fangbo Pan and Jinmei Pan, September 25, 2023
% jiang@bnu.edu.cn,panfb@mail.bnu.edu.cn
%% ********************************************************************* %%

clear
clc
close all

cTable = defaultSnowColormap;
TiffTags = struct('Compression',Tiff.Compression.LZW);

load('DEM_005_ASTER_AsiaTower.mat');% Auxiliary data: DEM
DEMdata=DEM_005_ASTER;


latLimits = [24,54];
lonLimits = [60,106];
res=0.005;
lon_005=[lonLimits(1)+res/2:res:lonLimits(2)];
lat_005=[latLimits(2)-res/2:-res:latLimits(1)];
refObj = georefcells(latLimits,lonLimits,[length(lat_005),length(lon_005)],'ColumnsStartFrom','north');

infolder='/home/Panfb/AWT_MODIS_FSC/CloudMasked_data/';
outfolder='/home/Panfb/AWT_MODIS_FSC/GapFilled_data/';

%Keys:
KeyWater=237;
KeyCloud=250;
KeyBoundary=255;
    
for iyear=1999:2022
    timebuffer=30; 
    % the period from July 1 to June 30 of the following year is taken as a year of snow cover
    if(iyear==1999)
        dates0=datenum(2000,2,26):1:datenum(2000,6,30)-1;
        dates=datenum(2000,2,26):1:datenum(2000,6,30)-1+timebuffer;
    elseif(iyear==2022)
        dates0=datenum(2022,7,1):1:datenum(2022,12,31);
        dates=datenum(2022,7,1)-timebuffer:1:datenum(2022,12,31);
    else
        dates0=datenum(iyear,7,1):1:datenum(iyear+1,6,30)-1;
        dates=datenum(iyear,7,1)-timebuffer:1:datenum(iyear+1,6,30)-1+timebuffer;
    end
    clear infiles oufiles
    ndates=length(dates);
    
    %Define the file to save
    data=zeros(size(DEMdata,1),size(DEMdata,2),length(dates),'single');
    for id=1:length(dates)
        [year,month,day]=datevec(dates(id));
        doy=dates(id)-datenum(year,1,0);
        yearstr=num2str(year);
        doystr=num2str(year*1000+doy); doystr=doystr(5:7);
        infiles{id}=[infolder,yearstr,'\MCDAGE_C6_',yearstr,'_',doystr,'.fSCA.color.tif'];
        infiles{id}=[infolder,'AWT_MODIS_FSC_',datestr(dates(id),'yyyymmdd'),'.tif'];
        try 
            disp(infiles{id})
            Snowtemp = geotiffread(infiles{id});
            Snowtemp(Snowtemp>250) = 250;
            data(:,:,id) = Snowtemp;
        catch
            data(:,:,id)=zeros(size(DEMdata,1),size(DEMdata,2))+KeyCloud;
        end
        clear Snowtemp
    end
      
    %% Step0:Record the the lake and nodata
    lakemask=sum(data==237,3)>10;
    for i=1:ndates
        tmp=data(:,:,i);
        tmp(tmp==KeyWater & lakemask==0)=KeyCloud;
        data(:,:,i)=tmp;
    end
        
    a=sum((data==0 | data==255),3);
    cr=max(a,[],'includenan');%,'omitnan');
    boundarymask=(a==cr);
    clear a cr
   %% Step1: Temporal interpolation
    for i=2:size(data,3)-1
        subdataIn=data(:,:,[i-1,i+1]);
        idx=find(subdataIn>100); subdataIn(idx)=NaN;
        subdataInmean=mean(subdataIn,3);
        subdataOut=data(:,:,i);
        idx=find(subdataOut==KeyCloud & isnan(subdataInmean)==0);
        if(length(idx)>0)
            subdataOut(idx)=subdataInmean(idx);
        end
        subdataOut(lakemask)=KeyWater;
        subdataOut(boundarymask)=KeyBoundary;
        data(:,:,i)=subdataOut;
    end   
    %% Step2: Spatio interpolation
    for i=1:ndates
         tmp=data(:,:,i);
         tmpS=nan(size(data,1),size(data,2),4);
         tmpS(2:end-1,2:end-1,1)=tmp(1:end-2,1:end-2);
         tmpS(2:end-1,2:end-1,2)=tmp(3:end,3:end);
         tmpS(2:end-1,2:end-1,3)=tmp(1:end-2,3:end);
         tmpS(2:end-1,2:end-1,4)=tmp(3:end,1:end-2);
         idx=find(tmpS>100);
         if(length(idx)>0); tmpS(idx)=NaN; end
         meanValue=mean(tmpS,3,'OmitNaN');
         ctSnow=sum(tmpS>0,3);
         idx=find(ctSnow>3 & tmp==KeyCloud);
         if(length(idx)>0)
             tmp(idx)=meanValue(idx);
         end
         ctLand=sum(tmpS==0,3);
         idx=find(ctLand>3 & tmp==KeyCloud);
         if(length(idx)>0)
             tmp(idx)=0;
         end
         data(:,:,i)=tmp;
    end
   %% Step3: Temporal interpolation
    x=1:ndates;
    for i=1:size(data,1)
        for j=1:size(data,2)     
            if(lakemask(i,j)==1 | boundarymask(i,j)==1)
                continue
            end                
            subdata=squeeze(data(i,j,:));
            idx=find(subdata==KeyCloud);
            if(length(idx)==0)
                continue
            end            
            subdataout=subdata;
            idx=find(subdata<=100);
            idxn=find(subdata==KeyCloud);
            if(length(idx)>2)
                subdataout(idxn)=interp1(x(idx),subdata(idx),x(idxn),'linear');
                for k=10:ndates-9
                    if(subdata(k)~=KeyCloud)
                        continue
                    end
                    tmp=subdata(k-9:k+9);
                    idx=find(tmp<=100);
                    if(length(idx)==0)
                        subdataout(k)=KeyCloud;
                    end
                end
            end
            idxLake=find(subdata==KeyWater);
            idxBoundary=find(subdata==KeyBoundary);
            if(length(idxLake)>0)
                subdataout(idxLake)=KeyWater;
            end
            if(length(idxBoundary)>0)
                subdataout(idxBoundary)=KeyBoundary;
            end            
            data(i,j,:)=subdataout;
        end
    end
   %% Step4: Spatio interpolation     
    for i=1:ndates
        if(dates(i)<dates0(1) | dates(i)>dates0(end))
            continue
        end
        
        tmp=squeeze(data(:,:,i));
        
        idxn=find(tmp==KeyCloud);
        if(length(idxn)>0)
            tmp(idxn)=NaN;
        else
            continue
        end
        tmp2=tmp;
        sizetmp=size(tmp);
        for k=1:length(idxn)
            [i0,j0]=ind2sub(sizetmp,idxn(k));
            ii=max(i0-5,1):min(i0+5,sizetmp(1));
            jj=max(j0-5,1):min(j0+5,sizetmp(2));
            blk=double(tmp(ii,jj));
            blk_dem = DEMdata(ii,jj);%dem
            idx=find(blk>100);
            if(length(idx)>0)
                blk(idx)=NaN;
            end
            if( sum(~isnan(blk),'all','omitnan')==0)
                continue
            else
                ix = abs(blk_dem-blk_dem(6,6))>100;
                blk(ix) = nan;
                for m = 1:11
                    for n = 1:11
                        dist(m,n) = sqrt((m-6).^2+(n-6).^2);
                    end
                end
                temp_dist = dist_dem./dist;
                temp_dist(6,6) = 0;
                w = temp_dist/sum(sum(temp_dist));
                idwvalue = sum(sum(blk.*w));
                tmp2(idxn(k))=idwvalue;
            end
        end            
        tmp2(lakemask)=KeyWater;
        tmp2(boundarymask)=KeyBoundary;
        idx=find(isnan(tmp2)==1);
        tmp2(idx)=KeyCloud;
        data(:,:,i)=tmp2;
    end
   %% Tiff file output       
        for i=1:ndates
            if(dates(i)<dates0(1) | dates(i)>dates0(end))
                continue
            end        
            filenameout=outfiles{i};
            disp(filenameout);
            dataout=data(:,:,i);
            geotiffwrite(filenameout,uint8(dataout),cTable,refObj,'TiffTags',TiffTags);
        end     
end