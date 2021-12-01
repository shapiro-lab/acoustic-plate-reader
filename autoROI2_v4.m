close all

showf5 = 1; % 0 or 1 to display ROIs on images
savedata = 1; % 0 or 1 to save processed data
xROI_size = 10;
zROI_size = 40;
zoffset = 80;
testfilter = 0;
rangeoffset = 10;
%%
sampROI = nan(Nf,2,total_n);
noiseROI_mean = sampROI;
noiseROI_std = sampROI;
nz = nan(2,total_n);
dz = nan(2,total_n);
%% manual correction of the ROI selection
xcorrection = zeros(1,total_n);
zcorrection = zeros(1,total_n);
% xcorrection(13) = -10;
% xcorrection(15) = -5;
% zcorrection(13) = 30;
% zcorrection(21:23) = 25;
% zcorrection(24) = 35;
%%
for testInd = 1:total_n
    Imt = Imi{1,2,testInd};
    filter_size1 = [40 5];

    threshold_ratio = 1;
    noise_stdratio = 1;
    auto_noise = 1;
    sample_depth = [3 8];
    noise_slice = [2 3];
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

    Imt_f1 = medfilt2(Imt,filter_size1);
    threshold1 = noiseROIt*threshold_ratio;
    Imt_f2 = (Imt_f1 > threshold1) .* filter_depth;
    if testfilter
        figure;
        imagesc(Xi, Zi, Imt_f2); axis image;
    end
    
    Xm = repmat(1:length(Xi),length(Zi),1).* Imt_f2;
    Xm2 = Imt_f2;
    % Xm(Xm==0) = nan;
    % bx1 = max(Xm,[],2) - min(Xm,[],2);
    bx2 = sum(Xm2,2);
    [~,ztop] = max(bx2);
    Ind_dZ = [ztop-zROI_size ztop+zROI_size];
    Ind_dZ = Ind_dZ + zoffset + zcorrection(testInd);
    ixcenter = round(sum(sum((Xm(Ind_dZ(1):Ind_dZ(2),:)))) / nnz(Xm(Ind_dZ(1):Ind_dZ(2),:)));
    Ind_dX = [ixcenter-xROI_size ixcenter+xROI_size]+ xcorrection(testInd);
    if showf5
        fig = figure;
        imagesc(Xi, Zit, 20*log10(abs(Imt(iZt,:))), [20 80]); axis image;colormap hot
        hold on;
        drawrectangle('Position',[Xi(Ind_dX(1)) Zi(Ind_dZ(1)) Xi(Ind_dX(2))-Xi(Ind_dX(1)) Zi(Ind_dZ(2))-Zi(Ind_dZ(1))],'EdgeColor','w','LineWidth',2)
        pause(1);
        % close(fig);
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
if savedata
    save([saveName '_data_' datestr(now,'yymmdd-hh-MM-ss')],'sampROI','sampCNR','noiseROI_mean','noiseROI_std','voltage','PlateCoordinate');
end
%close all