   %overlay = 1;
showf5 = 0;
xROI_size = 10;
zROI_size = 40;
zoffset = 80;
testfilter = 0;
rangeoffset = 10;

sampROI = nan(Nf,2,total_n);
sampCNR = sampROI;
noiseROI_mean = sampROI;
noiseROI_std = sampROI;
nz = nan(2,total_n);
dz = nan(2,total_n);
for testInd = 1
    %Ind = 59;    
    Imt = Im{1,2,testInd};
    
    %dImi(:,:,imgMode,Ind);
    filter_size1 = [20 5];
    filter_size2 = [5 5];
    filter_sigma = 2.5;
    filter_box = [9 1];
    %SizeThreshold = 300;

    threshold_ratio = 1;
    noise_stdratio = 1;
    auto_noise = 1;
    noiseSize = [2 2];
    sample_depth = [3 8];
    noise_slice = [2 3];
    sample_width = 10;
    dz(:,testInd) = [find(Zi>=sample_depth(1),1,'first') find(Zi>=sample_depth(2),1,'first')];
    nz(:,testInd) = [find(Zi>=noise_slice(1),1,'first') find(Zi>=noise_slice(2),1,'first')];
    iZt = dz(1,testInd):dz(2,testInd);
    iZd_noise = nz(1,testInd):nz(2,testInd);
    filter_depth = zeros(size(Imt));
    filter_depth(iZt,:) = 1;
    Zit = Zi(iZt);

    if auto_noise
        noiseROIt = mean(mean(Imt(iZd_noise,:))) + noise_stdratio * std2(Imt(iZd_noise,:));
    else
        noiseROIt = noiseROI(:,Ind);
    end

    Imt_f1 = medfilt2(Imt(iZt,:),filter_size1);
    threshold1 = noiseROIt*1;
    Imt_f2 = (Imt_f1 > threshold1);
    if testfilter
        figure;
        imagesc(Xi, Zit, Imt_f2); axis image;
    end
    
    Xm = repmat(1:length(Xi),length(Zit),1).* Imt_f2;
    bx1 = max(Xm,[],2) - min(Xm,[],2);
    bx2 = sum(Imt_f2,2);
    [~,ztop] = max(bx1);
    izbot = max(any(Imt_f2,2).*(1:length(Zit))');
    izbot_range = -ztop + izbot - rangeoffset;
    ixcenter = round(sum(sum((Xm(izbot - izbot_range:izbot,:)))) / nnz(Xm(izbot - izbot_range:izbot,:)));
    
    Ind_dX = [ixcenter-xROI_size ixcenter+xROI_size];
    
    
    Ind_dZ = [ztop-zROI_size ztop+zROI_size];
    Ind_dZ = iZt(Ind_dZ);
    Ind_dZ = Ind_dZ + zoffset;
    
    if showf5
        fig = figure;
        imagesc(Xi, Zit, 20*log10(abs(Imt(iZt,:))), [20 80]); axis image;colormap hot
        hold on;
        rectangle('Position',[Xi(Ind_dX(1)) Zi(Ind_dZ(1)) Xi(Ind_dX(2))-Xi(Ind_dX(1)) Zi(Ind_dZ(2))-Zi(Ind_dZ(1))],'EdgeColor','w','LineWidth',2)
        pause(1);
        %close(fig);
    end
    
    for j = 1:Nf
        for k = 1:2
            sampROI(j,k,testInd) = mean(mean(Imi{j,k,testInd}(Ind_dZ(1):Ind_dZ(2),Ind_dX(1):Ind_dX(2))));
            noiseROI_mean(j,k,testInd) = mean(mean(Imi{j,k,testInd}(iZd_noise,:)));
            noiseROI_std(j,k,testInd) =  std2(Imi{j,k,testInd}(iZd_noise,:));
        end
    end
end
sampCNR = 20 * log10(abs(sampROI - noiseROI_mean) ./ noiseROI_std);
clear i j k;