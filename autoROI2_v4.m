close all

showf5 = 1; % 0 or 1 to display ROIs on images
savedata = 1; % 0 or 1 to save processed data
xROI_size = 10;
zROI_size = 40;
zoffset = 100;
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
skip = zeros(1,total_n);
noise_slices = repmat([2 3],total_n,1);
noise_slices_backup = [8 9;1 2;7 8];
ntt = 100; % threshold for choosing noise depth slice (mV)

%fig_n = figure; hold on;]
% xcorrection([59]) = -5;
% xcorrection([20]) = -5;
% xcorrection([22]) = -5;
% xcorrection([23]) = -5;
% xcorrection([32]) = -7;
% xcorrection([47]) = -5;
% xcorrection([49]) = 5;
% xcorrection([55]) = 5;
% xcorrection([56]) = -5;
% xcorrection([58]) = -5;
% xcorrection([67]) = -5;
% xcorrection([70]) = -5;
% xcorrection([72]) = -5;
% xcorrection([91]) = -5;
% xcorrection([95]) = -9;
% xcorrection([96]) = 3;

% zcorrection(21:24) = 25;
% zcorrection(13) = 2;
% zcorrection(22) = -2;
%% quantify ROIs
for wellInd = 1:total_n
    if ~skip(wellInd)
        Imt = Imi{1,2,wellInd};
        filter_size1 = [40 5];

        threshold_ratio = 1;
        noise_stdratio = 1.5;
        auto_noise = 1;
        sample_depth = [2.5 8];
        noise_slice = noise_slices(wellInd,:);
        dz(:,wellInd) = [find(Zi>=sample_depth(1),1,'first') find(Zi>=sample_depth(2),1,'first')];
        nz(:,wellInd) = [find(Zi>=noise_slice(1),1,'first') find(Zi>=noise_slice(2),1,'first')];
        iZt = dz(1,wellInd):dz(2,wellInd);
        iZd_noise = nz(1,wellInd):nz(2,wellInd);
        filter_depth = zeros(size(Imt));
        filter_depth(iZt,:) = 1;
        Zit = Zi(iZt);
        noiseROIt = mean(mean(Imt(iZd_noise,:))) + noise_stdratio * std2(Imt(iZd_noise,:));
        i_ntt = 1;
        while noiseROIt > ntt && (i_ntt < length(noise_slices_backup(:,1)))
            noise_slice = noise_slices_backup(i_ntt,:);
            nz(:,wellInd) = [find(Zi>=noise_slice(1),1,'first') find(Zi>=noise_slice(2),1,'first')];
            iZd_noise = nz(1,wellInd):nz(2,wellInd);
            noiseROIt = mean(mean(Imt(iZd_noise,:))) + noise_stdratio * std2(Imt(iZd_noise,:));
            i_ntt = i_ntt + 1;
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
        Xm(Xm==0) = nan;
        bx1 = max(Xm,[],2) - min(Xm,[],2);
        bx2 = sum(Xm2,2);
        [~,ztop] = max(bx1);
        Ind_dZ = [ztop-zROI_size ztop+zROI_size];
        Ind_dZ = Ind_dZ + zoffset + zcorrection(wellInd);
        Xm(isnan(Xm)) = 0;
        ixcenter = round(sum(sum((Xm(Ind_dZ(1):Ind_dZ(2),:)))) / nnz(Xm(Ind_dZ(1):Ind_dZ(2),:)));
        if ~isnan(ixcenter)
            Ind_dX = [ixcenter-xROI_size ixcenter+xROI_size]+ xcorrection(wellInd);
            if showf5
                fig = figure;
                imagesc(Xi, Zi(1:iZt(end)+200), 20*log10(abs(Imt(1:iZt(end)+200,:))), [20 80]); axis image;colormap hot
                hold on;
                drawrectangle('Position',[Xi(1) Zi(iZd_noise(1)) Xi(end)-Xi(1) Zi(iZd_noise(end))-Zi(iZd_noise(1))],'EdgeColor','w','LineWidth',2);
                drawrectangle('Position',[Xi(Ind_dX(1)) Zi(Ind_dZ(1)) Xi(Ind_dX(2))-Xi(Ind_dX(1)) Zi(Ind_dZ(2))-Zi(Ind_dZ(1))],'EdgeColor','g','LineWidth',2);
                title(['Frame #' num2str(wellInd)]);
                pause(0.3);
                % close(fig);
            end

            for pressure = 1:Nf
                for imMode = 1:2
                    sampROI(pressure,imMode,wellInd) = mean(mean(Imi{pressure,imMode,wellInd}(Ind_dZ(1):Ind_dZ(2),Ind_dX(1):Ind_dX(2))));
                    noiseROI_mean(pressure,imMode,wellInd) = mean(mean(Imi{pressure,imMode,wellInd}(iZd_noise,:)));
                    noiseROI_std(pressure,imMode,wellInd) =  std2(Imi{pressure,imMode,wellInd}(iZd_noise,:));
                end
            end
        else
            ixcenter = round(length(Xi)/2);
            Ind_dX = [ixcenter-xROI_size ixcenter+xROI_size]+ xcorrection(wellInd);
            if showf5
                fig = figure;
                imagesc(Xi, Zi(1:iZt(end)), 20*log10(abs(Imt(1:iZt(end),:))), [20 80]); axis image;colormap hot
                hold on;
                drawrectangle('Position',[Xi(1) Zi(iZd_noise(1)) Xi(end)-Xi(1) Zi(iZd_noise(end))-Zi(iZd_noise(1))],'EdgeColor','w','LineWidth',2);
                drawrectangle('Position',[Xi(Ind_dX(1)) Zi(Ind_dZ(1)) Xi(Ind_dX(2))-Xi(Ind_dX(1)) Zi(Ind_dZ(2))-Zi(Ind_dZ(1))],'EdgeColor','g','LineWidth',2);
                title(['Warning: ROI not detected in frame #' num2str(wellInd)]);
                pause(0.3);
                % close(fig);
            end

            for pressure = 1:Nf
                for imMode = 1:2
                    sampROI(pressure,imMode,wellInd) = mean(mean(Imi{pressure,imMode,wellInd}(Ind_dZ(1):Ind_dZ(2),Ind_dX(1):Ind_dX(2))));
                    noiseROI_mean(pressure,imMode,wellInd) = mean(mean(Imi{pressure,imMode,wellInd}(iZd_noise,:)));
                    noiseROI_std(pressure,imMode,wellInd) =  std2(Imi{pressure,imMode,wellInd}(iZd_noise,:));
                end
            end
        end
    end
end
sampCNR = 20 * log10(abs(sampROI - noiseROI_mean) ./ noiseROI_std);

% save data
clear i pressure imMode;
if savedata
    save([saveName '_data_' datestr(now,'yymmdd-hh-MM-ss')],'P','saveName','sampROI','sampCNR','noiseROI_mean','noiseROI_std','voltage','PlateCoordinate');
end

%% plot ROI quants with microplateplot
% reshape ROI CNRs
% sampCNR new dimensions: well rows, well columns, frames, imaging modes
sampCNRs = permute(reshape(sampCNR, Nf, 2, PlateSize(2), PlateSize(1)), [4 3 1 2]);

% find max signal achieved by each sample at any voltage
maxs = squeeze(max(sampCNRs, [], 3));

% make and save microplate plot
mpplot = microplateplot(maxs(:,:,1));
colormap hot
colorbar
title('Max signal achieved at any voltage')
mpplot
savefig([saveName '_max-signal'])