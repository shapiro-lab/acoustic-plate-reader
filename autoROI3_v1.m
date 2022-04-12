close all
clear iZt

sample_depth = [1.9 9.1]; % display/sample dpeth range in mm
DisplayMode = 1; % How to display images: 0 = show nothing, 1 = always show images + ROIs, 2 = only show images + ROIs with inconfident ROI selection
testfilter = 0; % 0 or 1 to display the thresholding filter
savedata = 1; % 0 or 1 to save processed data
xROI_size = 10; % size of the ROI in x dimension    
zROI_size = 40; % size of the ROI in z dimension 
filter_size1 = [40 5]; % size of the median filter in [z,x] dimension
noise_stdratio = 2; % ratio of the noise STD added to the noise mean to decide noise floor 
TemplateMode = 3; % ROI template to use: 1 = filled wells, 2/3 = empty wells with thinner (2) or thicker (3) well wall, use 1 as default and 2/3 if well is not filled
if ~exist('WellTemplates','var') % load the ROI template if needed, please put the .mat file in accessible folers
    load('WellTemplates.mat');
end
WellTemplate = WellTemplates{TemplateMode}; % Apply ROI template as selected.
CST = ConfidenceScoreThreshold(TemplateMode); % Apply corresponding confidence score threshold for reminder to adjust ROIs (abitrary defined thresholds)
%%
sampROI = nan(Nf,2,total_n); % initialize for sample ROI
noiseROI_mean = sampROI; % initialize for noise ROI
noiseROI_std = sampROI; % initialize for noise ROI STD
nz = nan(2,total_n); % initialize for noise slice positions
ConfidenceScore = nan(1,total_n); % initialize for ROI prediction confidence score 
ROI_Centers = nan(total_n,2); % % initialize for locations of ROI centers
%% manual correction of the ROI selection
xcorrection = zeros(1,total_n); % initialize for ROI correction in x dimension (correction in the unit of voxel count, ~0.1 mm/count)
zcorrection = 20*ones(1,total_n);% initialize for ROI correction in z dimension (correction in the unit of voxel count, ~0.0124 mm/count)  
skip = zeros(1,total_n); % ROI correction in well indices to skip processing
noise_slices = repmat([2 3],total_n,1); % default noise depth (in mm) for noise ROI selection
noise_slices_backup = [8 9;1 2;7 8]; % backup noise depth (in mm) for noise ROI selection in case of abnormal noise level
ntt = 100; % threshold of noise level (mV) for abnormal high noise level / changing noise depth slices

xcorrection([2]) = 2;
xcorrection([13]) = 5;
xcorrection([22]) = 5;
xcorrection([31]) = 5;

% xcorrection([32]) = -7;
% xcorrection([47]) = -5;
% xcorrection([12]) = -5;
% xcorrection([18]) = 5;
% xcorrection([34]) = -5;
% xcorrection([38]) = 5;
% xcorrection([60]) = -5;
% xcorrection([27]) = 5;
% xcorrection([53]) = -5;
% xcorrection([55]) = -5;
% xcorrection([78]) = 5;

% zcorrection(21:24) = 25;
% zcorrection(7) = -40;
% zcorrection(22) = -2;

%% quantify ROIs
for wellInd = 1:total_n
    if ~skip(wellInd)
        Imt = Imi{1,2,wellInd}; % Take the Bmode at the first voltage to detect ROI
        noise_slice = noise_slices(wellInd,:); % Fetch noise depth slice to select noise ROI
        if ~exist('iZt','var') % Select display/sample depth if not done yet
            iZt = find(Zi>=sample_depth(1),1,'first'):find(Zi>=sample_depth(2),1,'first'); 
        end
        Imt = Imt(iZt,:); % Crop the images to defined sample depth for computational efficiency
        Zit = Zi(iZt); % Define new z axis after croping 
        nz(:,wellInd) = [find(Zit>=noise_slice(1),1,'first') find(Zit>=noise_slice(2),1,'first')]; %Find start/end indices for corresponding noise ROI
        iZd_noise = nz(1,wellInd):nz(2,wellInd); % Convert to index array
        noiseROIt = mean(mean(Imt(iZd_noise,:))) + noise_stdratio * std2(Imt(iZd_noise,:)); % Calculate the noise level for thresholding 
        i_ntt = 1; % Initialize counting for calculating based on backup noise depth slices 
        while noiseROIt > ntt && (i_ntt < length(noise_slices_backup(:,1))) % Change noise depth slices until below the threshold / running out slices
            noise_slice = noise_slices_backup(i_ntt,:); % Fetch new slice
            nz(:,wellInd) = [find(Zit>=noise_slice(1),1,'first') find(Zit>=noise_slice(2),1,'first')]; % Convert to indices
            iZd_noise = nz(1,wellInd):nz(2,wellInd);% Convert to indices
            noiseROIt = min(noiseROIt,mean(mean(Imt(iZd_noise,:))) + noise_stdratio * std2(Imt(iZd_noise,:)));% Take lower one between the old and new noise level 
            i_ntt = i_ntt + 1; % Iterate 
        end
        
        Imt_f1 = medfilt2(Imt,filter_size1); % Apply median filter to remove pepper&salt noise
        Imt_f2 = (Imt_f1 > noiseROIt); % Apply thresholding based on the noise level
        if testfilter % Debugging mode, display applied filters
            figure;
            imagesc(Xi, Zit, Imt_f2); axis image;
        end
        c = xcorr2(WellTemplate,double(Imt_f2)); % Calculate 2D cross correlation between images and centered ROI template 
        [zpeak,xpeak] = find(c==max(c(:)),1); % Find the max correlation location
        ConfidenceScore(wellInd) = c(zpeak,xpeak); % Use the correlation at max location as confidence score
        zoffset = zpeak - size(Imt_f2,1) + zcorrection(wellInd); % Convert max correlation location to index offset in z dimension
        xoffset = xpeak - size(Imt_f2,2) - xcorrection(wellInd); % Convert max correlation location to index offset in z dimension
        ROI_Centers(wellInd,:) = TemplateCenter - [zoffset xoffset]; % Apply the offset to define the ROI center
        Ind_dZ = [ROI_Centers(wellInd,1) - zROI_size ROI_Centers(wellInd,1) + zROI_size]; % Extend from center to get start/end indices in z dimension 
        Ind_dX = [ROI_Centers(wellInd,2) - xROI_size ROI_Centers(wellInd,2) + xROI_size]; % Extend from center to get start/end indices in x dimension
        


        if DisplayMode == 1 % Always-on display mode, show images and predicted ROIs
            fig = figure;
            imagesc(Xi, Zit, 20*log10(abs(Imt)), [20 80]); axis image;colormap hot
            hold on;
            drawrectangle('Position',[Xi(1) Zit(iZd_noise(1)) Xi(end)-Xi(1) Zit(iZd_noise(end))-Zit(iZd_noise(1))],'EdgeColor','w','LineWidth',2);
            drawrectangle('Position',[Xi(Ind_dX(1)) Zit(Ind_dZ(1)) Xi(Ind_dX(2))-Xi(Ind_dX(1)) Zit(Ind_dZ(2))-Zit(Ind_dZ(1))],'EdgeColor','g','LineWidth',2);
            title(['Frame #' num2str(wellInd)]);
            pause(0.3);
            % close(fig);
        elseif DisplayMode == 2 && (ConfidenceScore(wellInd) < CST)% Only below confidence score threshold images would be displayed 
            fig = figure;
            imagesc(Xi, Zit, 20*log10(abs(Imt)), [20 80]); axis image;colormap hot
            hold on;
            drawrectangle('Position',[Xi(1) Zit(iZd_noise(1)) Xi(end)-Xi(1) Zit(iZd_noise(end))-Zit(iZd_noise(1))],'EdgeColor','w','LineWidth',2);
            drawrectangle('Position',[Xi(Ind_dX(1)) Zit(Ind_dZ(1)) Xi(Ind_dX(2))-Xi(Ind_dX(1)) Zit(Ind_dZ(2))-Zit(Ind_dZ(1))],'EdgeColor','g','LineWidth',2);
            title(['Frame #' num2str(wellInd)]);
            pause(0.3);
        end

        for pressure = 1:Nf % Calulate sample ROI, noise ROI, and noise STD using the predicted ROI for all the voltages
            for imMode = 1:2
                Imt = Imi{pressure,imMode,wellInd}(iZt,:);
                sampROI(pressure,imMode,wellInd) = mean(mean(Imt(Ind_dZ(1):Ind_dZ(2),Ind_dX(1):Ind_dX(2))));
                noiseROI_mean(pressure,imMode,wellInd) = mean(mean(Imt(iZd_noise,:)));
                noiseROI_std(pressure,imMode,wellInd) =  std2(Imt(iZd_noise,:));
            end
        end
    end
end
sampCNR = 20 * log10(abs(sampROI - noiseROI_mean) ./ noiseROI_std); % Calculate sample CNR

% save data
clear pressure imMode;
if savedata
    save([saveName '_data_' datestr(now,'yymmdd-hh-MM-ss')],'P','saveName','sampROI','sampCNR','noiseROI_mean','noiseROI_std','voltage','PlateCoordinate','ROI_Centers','ConfidenceScore','Nf');
end

%% plot ROI quants with microplateplot
% reshape ROI CNRs
% sampCNR new dimensions: well rows, well columns, frames, imaging modes
sampCNRs = permute(reshape(sampCNR, Nf, 2, PlateSize(2), PlateSize(1)), [4 3 1 2]);
ConfidenceScore = permute(reshape(ConfidenceScore, PlateSize(2), PlateSize(1)), [1 2]);

% find max signal achieved by each sample at any voltage
maxs = squeeze(max(sampCNRs, [], 3));

% make and save microplate plot
figure;
mpplot = microplateplot(maxs(:,:,1));
colormap hot
colorbar
title('Max signal achieved at any voltage')
mpplot
savefig([saveName '_max-signal'])