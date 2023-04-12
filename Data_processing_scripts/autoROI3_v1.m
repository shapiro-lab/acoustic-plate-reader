% clear all
% close all
% 
% pathName = 'G:\.shortcut-targets-by-id\0B24ONICaZ0z9djczVE1ZR3BnWU0\Shapiro Lab Information\Data\Rob\96-well_plate_scans\GvpA-B-mutants\A-lib-2\A-lib-K22R-A2\A-lib-K22R-A2_P3_4_stable-37C_P_R1_C3';

% PlateProc1_v2
% 
close all
clear iZt ixZtemp Im
%%
sample_depth = [1 10]; % display depth range in mm
ROI_depth = [3 8]; % sample depth range in mm
Bmode_color_lims = [0 60];
xAM_color_lims = [0 70];

%Mohamed 12/27: DisplayMode is now obsolete, only thing that shows is layout and 96 well plots
DisplayMode = 1; % How to display images: 0 = show nothing, 1 = always show images + ROIs, 2 = only show images + ROIs with low-confidence ROI selection

testfilter = 0; % 0 or 1 to display the thresholding filter
xROI_size = 10; % size of the ROI in x dimension    
zROI_size = 40; % size of the ROI in z dimension 
filter_size1 = [40 5]; % size of the median filter in [z,x] dimension
noise_stdratio = 3; % ratio of the noise STD added to the noise mean to decide noise floor 
noise_slices = repmat([8 9],total_n,1); % default noise depth (in mm) for noise ROI selection
noise_slices_backup = [2 3;9 10;10 11]; % backup noise depth (in mm) for noise ROI selection in case of abnormal noise level
noiseThresh = 100; % threshold of noise level (mV) for abnormal high noise level / changing noise depth slices
TemplateMode = 3; % ROI template to use: 1 = filled wells, 2/3 = empty wells with thinner (2) or thicker (3) well wall, use 1 as default and 2/3 if well is not full
if ~exist('WellTemplates','var') % load the ROI template if needed, please put the .mat file in accessible folders
    load('WellTemplates.mat');
end
WellTemplate = WellTemplates{TemplateMode}; % Apply ROI template as selected.
CST = ConfidenceScoreThreshold(TemplateMode); % Apply corresponding confidence score threshold for reminder to adjust ROIs (abitrary defined thresholds)
VmaxIX = find(P.Vseq==max(P.Vseq)); % find indices of all occurences of max voltage

%% the following arrays are used for callback purposes
Bmode_figs = gobjects(1,total_n); % initialize array to hold Bmode figure objects
xAM_figs = gobjects(1,total_n); % initialize array to hold xAM figure objects
ZixROI_array = nan(total_n,2); % initialize array to hold ZiX ROI values
XixROI_array = nan(total_n,2); % initialize array to hold XiX ROI values
Xi_cell = cell(total_n,1); %cell to hold Xi doubles
ZixTemp_cell = cell(total_n,1); %cell to hold ZiXTemp doubles
maxs_plots = gobjects(1,2); % initialize array to hold maxs_XAM and Bmode figures

%% other arrays that need to be initialized
noiseZ = nan(2,total_n); % initialize array for noise slice bounds
confScore = nan(1,total_n); % initialize array for ROI prediction confidence score 
ROI_Centers = nan(total_n,2); % initialize array for locations of ROI centers
quants.sampROI_means = nan(Nf,2,total_n); % initialize array for sample means
quants.noiseROI_means = nan(Nf,2,total_n); % initialize array for noise means
quants.noiseROI_stds = nan(Nf,2,total_n); % initialize array for noise ROI STDs
Ixz_cell = cell(Nf,1); %initialize empty cell to hold ixZtemp from each well

%% manual correction of ROI selection
xcorrection = zeros(1,total_n); % initialize array for ROI correction in x dimension (correction in the unit of voxel count, ~0.1 mm/count)
zcorrection = zeros(1,total_n); % initialize array for ROI correction in z dimension (correction in the unit of voxel count, ~0.0124 mm/count)  
skip = zeros(1,total_n); % ROI correction in well indices to skip processing
zcorrection = zcorrection + 20;
%% quantify ROIs
h = waitbar(0,'Quantifying ROIs');
for wellIx = 1:total_n
    if ~skip(wellIx)
        ImTemp = Imi{1,2,wellIx}; % Fetch Bmode image at first voltage to use for ROI detection
        noise_slice = noise_slices(wellIx,:); % Fetch noise depth slice to select noise ROI
        if ~exist('ixZtemp','var') % Select display depth if not done yet
            ixZtemp = find(Zi>=sample_depth(1),1,'first'):find(Zi>=sample_depth(2),1,'first'); 
        end
        ImTemp = ImTemp(ixZtemp,:); % Crop the images to defined sample depth to decrease data size
        Ixz_cell(wellIx) = {ixZtemp}; % store ImTemp in cell for callback purposes
        ZixTemp = Zi(ixZtemp); % Define new z-axis after cropping
        if ~exist('ixZtemp_ROI','var') % Select display depth if not done yet
            ixZtemp_ROI = find(ZixTemp>=ROI_depth(1),1,'first'):find(ZixTemp>=ROI_depth(2),1,'first'); 
        end
        noiseZ(:,wellIx) = [find(ZixTemp>=noise_slice(1),1,'first') find(ZixTemp>=noise_slice(2),1,'first')]; % Find start/end indices for corresponding noise ROI
        iZd_noise = noiseZ(1,wellIx):noiseZ(2,wellIx); % Make vector of noise Z indices
        noiseROItemp = mean(mean(ImTemp(iZd_noise,:))) + noise_stdratio * std2(ImTemp(iZd_noise,:)); % Calculate the noise level for thresholding 
        i_noiseThresh = 1; % Initialize counting for calculating based on backup noise depth slices 
        while noiseROItemp > noiseThresh && (i_noiseThresh < length(noise_slices_backup(:,1))) % Change noise depth slices until below the threshold or run out of slices
            noise_slice = noise_slices_backup(i_noiseThresh,:); % Fetch new slice
            noiseZ(:,wellIx) = [find(ZixTemp>=noise_slice(1),1,'first') find(ZixTemp>=noise_slice(2),1,'first')]; % Make array of noise Z bounds
            iZd_noise = noiseZ(1,wellIx):noiseZ(2,wellIx); % Update vector of noise Z indices
            noiseROItemp = min(noiseROItemp,mean(mean(ImTemp(iZd_noise,:))) + noise_stdratio * std2(ImTemp(iZd_noise,:))); % Take lower one between the old and new noise level 
            i_noiseThresh = i_noiseThresh + 1; % Iterate 
        end
        noiseROItemp = 60;
        depth_filter = zeros(size(ImTemp));
        depth_filter(ixZtemp_ROI,:) = 1;
        ImTemp_filt1 = medfilt2(ImTemp,filter_size1); % Apply median filter to remove salt&pepper noise
        ImTemp_filt2 = (ImTemp_filt1 > noiseROItemp) .* depth_filter; % Apply thresholding based on the noise level
        if testfilter % Debugging mode, display applied filters
            figure;
            imagesc(Xi, ZixTemp, ImTemp_filt2); axis image;
        end
        cor = xcorr2(WellTemplate,double(ImTemp_filt2)); % Calculate 2D cross-correlation between images and centered ROI template 
        [zpeak,xpeak] = find(cor==max(cor(:)),1); % Find the max correlation location
        confScore(wellIx) = cor(zpeak,xpeak); % Use the correlation at max location as confidence score
        zoffset = zpeak - size(ImTemp_filt2,1) + zcorrection(wellIx); % Convert max correlation location to index offset in z dimension and apply manual z correction
        xoffset = xpeak - size(ImTemp_filt2,2) - xcorrection(wellIx); % Convert max correlation location to index offset in x dimension and apply manual x correction
        ROI_Centers(wellIx,:) = TemplateCenter - [zoffset xoffset]; % Apply the offset to define the ROI center
        ZixROI = [ROI_Centers(wellIx,1) - zROI_size ROI_Centers(wellIx,1) + zROI_size]; % Extend from center to get start/end indices in z dimension 
        XixROI = [ROI_Centers(wellIx,2) - xROI_size ROI_Centers(wellIx,2) + xROI_size]; % Extend from center to get start/end indices in x dimension
        if ZixROI(1) <= 0
            ZixROI(1) = 1;
            ZixROI(2) = 1+2*zROI_size;
            disp(['Warning: well#' num2str(wellIx) ' might be off.'])
        elseif ZixROI(2) > length(ZixTemp)
            ZixROI(2) = length(ZixTemp);
            ZixROI(1) = ZixROI(2)-2*zROI_size;
            disp(['Warning: well#' num2str(wellIx) ' might be off.'])
        end
        if XixROI(1) <= 0
            XixROI(1) = 1;
            XixROI(2) = 1+2*xROI_size;
            disp(['Warning: well#' num2str(wellIx) ' might be off.'])
        elseif XixROI(2) > length(Xi)
            XixROI(2) = length(Xi);
            XixROI(1) = XixROI(2)-2*xROI_size;
            disp(['Warning: well#' num2str(wellIx) ' might be off.'])
        end
        %store data needed for callbacks
        ZixROI_array(wellIx, :) = ZixROI; %store rectangle coords in array
        XixROI_array(wellIx, :) = XixROI; %store rectangle coords in array
        
        %make the figures
        if DisplayMode == 1 % Always-on display mode, show images and predicted ROIs
            %create Bmode figure
            Bmode_fig = figure();
            imagesc(Xi, ZixTemp, 20*log10(abs(ImTemp)), Bmode_color_lims); axis image; colormap bone;
            hold on;
            noiseROI = drawrectangle('Position',[Xi(1) ZixTemp(iZd_noise(1)) Xi(end)-Xi(1) ZixTemp(iZd_noise(end))-ZixTemp(iZd_noise(1))], 'EdgeColor','w', 'LineWidth',2, 'Tag',['noise_' num2str(wellIx)]); % Draw noise ROI
            sampROI = drawrectangle('Position',[Xi(XixROI(1)) ZixTemp(ZixROI(1)) Xi(XixROI(2))-Xi(XixROI(1)) ZixTemp(ZixROI(2))-ZixTemp(ZixROI(1))], 'EdgeColor','g', 'LineWidth',2, 'Tag',['sampl_' num2str(wellIx)]); % Draw sample ROI
            addlistener(noiseROI,'ROIMoved',@updateROI);
            addlistener(sampROI,'ROIMoved',@updateROI);
            title(['Bmode Frame #' num2str(wellIx)]);
            Bmode_figs(wellIx) = Bmode_fig; % store figure object in array
            set(Bmode_fig, 'visible', 'off'); % hide figure
%             pause(0.3);
           
            %create xAM figure
            xAM_fig = figure();
            ImTemp_xAM = Imi{VmaxIX(1),1,wellIx}; % Fetch xAM image at first occurence of max voltage
%             ImTemp_xAM = Imi{VmaxIX(1),1,wellIx}-Imi{VmaxIX(2),1,wellIx}; % Generate pre-post-collapse difference image at max voltage
            ImTemp_xAM = ImTemp_xAM(ixZtemp,:);
            imagesc(Xi, ZixTemp, 20*log10(abs(ImTemp_xAM)), xAM_color_lims); axis image; colormap hot;
            hold on;
            noiseROI = drawrectangle('Position',[Xi(1) ZixTemp(iZd_noise(1)) Xi(end)-Xi(1) ZixTemp(iZd_noise(end))-ZixTemp(iZd_noise(1))], 'EdgeColor','w', 'LineWidth',2, 'Tag',['noise_' num2str(wellIx)]); % Draw noise ROI
            sampROI = drawrectangle('Position',[Xi(XixROI(1)) ZixTemp(ZixROI(1)) Xi(XixROI(2))-Xi(XixROI(1)) ZixTemp(ZixROI(2))-ZixTemp(ZixROI(1))], 'EdgeColor','g', 'LineWidth',2, 'Tag',['sampl_' num2str(wellIx)]); % Draw sample ROI
            addlistener(noiseROI,'ROIMoved',@updateROI);
            addlistener(sampROI,'ROIMoved',@updateROI);
            title(['xAM Frame #' num2str(wellIx)]);
            xAM_figs(wellIx) = xAM_fig; % store figure object in array
            set(xAM_fig, 'visible', 'off'); % hide figure
            Xi_cell(wellIx) = {Xi}; %store Xi in cell
            ZixTemp_cell(wellIx) = {ZixTemp}; % store ZiXTemp in cell
%             pause(0.3);
            
        elseif DisplayMode == 2 && (confScore(wellIx) < CST) % Only display images below confidence score threshold
            Bmode_fig = figure(wellIx);
            imagesc(Xi, ZixTemp, 20*log10(abs(ImTemp)), Bmode_color_lims); axis image; colormap bone;
            hold on;
            noiseROI = drawrectangle('Position',[Xi(1) ZixTemp(iZd_noise(1)) Xi(end)-Xi(1) ZixTemp(iZd_noise(end))-ZixTemp(iZd_noise(1))], 'EdgeColor','w', 'LineWidth',2, 'Tag',['noise_' num2str(wellIx)]); % Draw noise ROI
            sampROI = drawrectangle('Position',[Xi(XixROI(1)) ZixTemp(ZixROI(1)) Xi(XixROI(2))-Xi(XixROI(1)) ZixTemp(ZixROI(2))-ZixTemp(ZixROI(1))], 'EdgeColor','g', 'LineWidth',2, 'Tag',['sampl_' num2str(wellIx)]); % Draw sample ROI
            addlistener(noiseROI,'ROIMoved',@updateROI);
            addlistener(sampROI,'ROIMoved',@updateROI);
            title(['Frame #' num2str(wellIx)]);
            Bmode_figs(wellIx) = Bmode_fig; % store figure object in array
            set(Bmode_fig, 'visible', 'off'); % hide figure
            %pause(0.3);
            
            %create xAM figure
            xAM_fig = figure();
            ImTemp_xAM = Imi{Nf-1,1,wellIx}; %Fetch xAM image, currently using first voltage image
            ImTemp_xAM = ImTemp_xAM(ixZtemp,:);
            imagesc(Xi, ZixTemp, 20*log10(abs(ImTemp_xAM)), xAM_color_lims); axis image; colormap hot;
            hold on;
            noiseROI = drawrectangle('Position',[Xi(1) ZixTemp(iZd_noise(1)) Xi(end)-Xi(1) ZixTemp(iZd_noise(end))-ZixTemp(iZd_noise(1))], 'EdgeColor','w', 'LineWidth',2, 'Tag',['noise_' num2str(wellIx)]); % Draw noise ROI
            sampROI = drawrectangle('Position',[Xi(XixROI(1)) ZixTemp(ZixROI(1)) Xi(XixROI(2))-Xi(XixROI(1)) ZixTemp(ZixROI(2))-ZixTemp(ZixROI(1))], 'EdgeColor','g', 'LineWidth',2, 'Tag',['sampl_' num2str(wellIx)]); % Draw sample ROI
            addlistener(noiseROI,'ROIMoved',@updateROI);
            addlistener(sampROI,'ROIMoved',@updateROI);
            title(['xAM Frame #' num2str(wellIx)]);
            xAM_figs(wellIx) = xAM_fig; % store figure object in array
            set(xAM_fig, 'visible', 'off'); % hide figure
            Xi_cell(wellIx) = {Xi}; %store Xi in cell
            ZixTemp_cell(wellIx) = {ZixTemp}; % store ZiXTemp in cell
            %pause(0.3);
        end

        for frame = 1:Nf % Calulate sample mean, noise mean, and noise STD using the predicted ROI for all frames
            for imMode = 1:2
                ImTemp = Imi{frame,imMode,wellIx}(ixZtemp,:); % get image for quantification
%                 ImTemp = Imi{frame,imMode,wellIx}(ixZtemp,:) - Imi{frame+13,imMode,wellIx}(ixZtemp,:); % get pre-post-collapse subtraction image for quantification
                quants.sampROI_means(frame,imMode,wellIx) = mean(mean(ImTemp(ZixROI(1):ZixROI(2),XixROI(1):XixROI(2))));
                quants.noiseROI_means(frame,imMode,wellIx) = mean(mean(ImTemp(iZd_noise,:)));
                quants.noiseROI_stds(frame,imMode,wellIx) = std2(ImTemp(iZd_noise,:));
            end
        end
    end
    waitbar(wellIx/total_n)
end
close(h);

quantify(quants, VmaxIX);

%% plot ROI quants with microplateplot
% make and save microplate plots
make_plots(quants,PlateSize,Nf,confScore,saveName)

sampSBR_reshaped = permute(reshape(quants.sampSBR, [], 2, PlateSize(2), PlateSize(1)), [4 3 1 2]);
sampSBR_diff_reshaped = permute(reshape(quants.sampSBR_diff, [], 2, PlateSize(2), PlateSize(1)), [4 3 1 2]);
means = squeeze(mean(sampSBR_diff_reshaped,1));

newcolors = {'#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#D95319', '#7f7f7f', '#bcbd22', '#17becf'};
% Plot best Round 2 mutants
figure;
% plot(squeeze(quants.AM_Bmode_ratio_diff(:,1,[23 9])));
% plot(squeeze(normalize(quants.AM_Bmode_ratio_diff(:,1,[14 9]),1,'range')));
% plot(squeeze(normalize(quants.sampSBR_diff(:,1,[5:9]),1,'range')));
% plot(squeeze(normalize(means([5:10],:,1)',1,'range')));
% plot(squeeze(means([14 18],:,1)'));
% plot(normalize(squeeze(sampSBR_diff_reshaped(1,[1 2 3],:,1))'));
plot(squeeze(sampSBR_reshaped(1:4,4,:,1))');
title('Best Round 2 GvpA Mutants')
xlabel('Voltage')
ylabel('xAM diff SBR')
% ylim([0 5])
colororder(newcolors)
legend

Bmode_layout = layoutfigures(Bmode_figs,PlateSize(1),PlateSize(2), 'Bmode layout', 'bone', Bmode_color_lims); %create Bmode figure layout
xAM_layout = layoutfigures(xAM_figs,PlateSize(1),PlateSize(2), 'xAM layout', 'hot', xAM_color_lims); %create xAM figure layout

%% function definitions
function quantify(quants, VmaxIX)
    % Function to produce quantifications of ROIs
    
    % calculate SBRs
    quants.sampSBR = quants.sampROI_means ./ quants.noiseROI_means; % Calculate sample SBR
    quants.sampSBR_dB = 20 * log10(quants.sampSBR); % Calculate sample SBR in dB
%     quants.sampSBR_diff = quants.sampSBR(2:VmaxIX(1),:,:) - quants.sampSBR(VmaxIX(1)+1:VmaxIX(2),:,:); % Calculate pre-/post-collapse difference SBR
%     quants.sampSBR_diff = quants.sampSBR(1:VmaxIX(1),:,:) - quants.sampSBR(VmaxIX(1)+1:VmaxIX(2),:,:); % Calculate pre-/post-collapse difference SBR
    quants.sampSBR_diff = quants.sampSBR(VmaxIX(1),:,:) - quants.sampSBR(VmaxIX(end),:,:); % Calculate pre-/post-collapse difference SBR
%     quants.sampSBR_diff = quants.sampSBR(1,:,:) - quants.sampSBR(2,:,:); % Calculate pre-/post-collapse difference SBR
    quants.sampSBR_diff_dB = 20 * log10(abs(quants.sampSBR_diff)); % Calculate pre-/post-collapse difference SBR in dB
    % plot(squeeze(quants.sampSBR(1:13,2,[9 21 33 45])))
    %plot(squeeze(quants.sampSBR(1:13,2,[10 22 34 46])))
    % calculate CNRs
    quants.sampCNR = (quants.sampROI_means - quants.noiseROI_means) ./ quants.noiseROI_stds; % Calculate sample CNR
    quants.sampCNR_dB = 20 * log10(abs(quants.sampCNR)); % Calculate sample CNR in dB
%     quants.sampCNR_diff = quants.sampCNR(2:VmaxIX(1),:,:) - quants.sampCNR(VmaxIX(1)+1:VmaxIX(2),:,:); % Calculate pre-/post-collapse difference CNR
%     quants.sampCNR_diff = quants.sampCNR(1:VmaxIX(1),:,:) - quants.sampCNR(VmaxIX(1)+1:VmaxIX(2),:,:); % Calculate pre-/post-collapse difference CNR
    quants.sampCNR_diff = quants.sampCNR(VmaxIX(1),:,:) - quants.sampCNR(VmaxIX(end),:,:); % Calculate pre-/post-collapse difference CNR
%     quants.sampCNR_diff = quants.sampCNR(1,:,:) - quants.sampCNR(2,:,:); % Calculate pre-/post-collapse difference CNR
    quants.sampCNR_diff_dB = 20 * log10(abs(quants.sampCNR_diff)); % Calculate pre-/post-collapse difference CNR in dB
    
    % calculate xAM/Bmode ratios
%     quants.AM_Bmode_ratio_diff = (quants.sampROI_means(1:VmaxIX(1),1,:) - quants.sampROI_means(VmaxIX(1)+1:VmaxIX(2),1,:)) ./ (quants.sampROI_means(1:VmaxIX(1),2,:) - quants.sampROI_means(VmaxIX(1)+1:VmaxIX(2),2,:)); % xAM/Bmode
    quants.AM_Bmode_ratio_diff = (quants.sampROI_means(VmaxIX(1),1,:) - quants.sampROI_means(VmaxIX(end),1,:)) ./ (quants.sampROI_means(VmaxIX(1),2,:) - quants.sampROI_means(VmaxIX(end),2,:)); % xAM/Bmode
%     quants.AM_Bmode_ratio_diff = (quants.sampROI_means(1:VmaxIX(1),1,:) - quants.sampROI_means(VmaxIX(1)+1:VmaxIX(2),1,:)) ./ (quants.sampROI_means(1,2,:) - quants.sampROI_means(VmaxIX(2),2,:)); % xAM/Bmode
%     quants.AM_Bmode_ratio = (quants.sampROI_means(:,1,:) - quants.noiseROI_means(:,1,:)) ./ (quants.sampROI_means(:,2,:) - quants.noiseROI_means(:,2,:)); % xAM/Bmode
    quants.AM_Bmode_ratio_dB = 20 * log10(abs(quants.AM_Bmode_ratio_diff)); % xAM/Bmode, dB scale
    
    assignin('base','quants',quants);
    
    % save updated quants
    evalin('base', "save([saveName '_data_' datestr(now,'yymmdd-hh-MM-ss')],'P','saveName','quants','voltage','PlateCoordinate','PlateSize','ROI_Centers','confScore','Nf');")
end

function make_plots(quants,PlateSize,Nf,confScore,saveName)
    % This function makes the maxs_AM and maxs_AM_Bmode_ratio figures. If they
    % exist, it deletes the previous maxs figures and recomputes the values
    
    %reshape ROI CNRs
    % sampCNR new dimensions: well rows, well columns, frames, imaging modes
    sampCNRs = permute(reshape(quants.sampSBR_diff, [], 2, PlateSize(2), PlateSize(1)), [4 3 1 2]);
    AM_Bmode_ratio_diffs = permute(reshape(quants.AM_Bmode_ratio_diff, length(quants.AM_Bmode_ratio_diff(:,1,1)), 1, PlateSize(2), PlateSize(1)), [4 3 1 2]);
    confScore = permute(reshape(confScore, PlateSize(2), PlateSize(1)), [1 2]);
    
    % find max signal achieved by each sample at any voltage
    maxs_AM = squeeze(max(sampCNRs, [], 3));
    maxs_AM_Bmode_ratio = squeeze(max(AM_Bmode_ratio_diffs, [], 3));
    
    %delete previous plots
    open_figs = findobj('type', 'figure');
    for fig_ind = 1:length(open_figs)
        fig = open_figs(fig_ind);
        if strcmp(fig.CurrentAxes.Tag, 'max-xAM') || strcmp(fig.CurrentAxes.Tag, 'max-xAM:Bmode') %this will delete the previous maxs figures if they exist, we remake them later
            delete(fig);
        end
    end
    
    % make and save microplate plots
    plot_xAM(maxs_AM, saveName);
    plot_AM_Bmode_ratio(maxs_AM_Bmode_ratio, saveName);
end

function plot_xAM(maxs_AM, saveName)
    figure;
    xAM_mpplot = microplateplot(maxs_AM(:,:,1));
    colormap hot
    colorbar
    title('Max xAM signal achieved at any voltage')
    set(xAM_mpplot, 'tag', 'max-xAM');
    xAM_mpplot;
    savefig([saveName '_max-xAM'])
end

function plot_AM_Bmode_ratio(maxs_AM_Bmode_ratio, saveName)
    figure;
    AM_Bmode_mpplot = microplateplot(maxs_AM_Bmode_ratio(:,:));
    colormap parula
    colorbar
    title('Max xAM:Bmode ratio signal achieved at any voltage')
    set(AM_Bmode_mpplot, 'tag', 'max-xAM:Bmode');
    AM_Bmode_mpplot;
    savefig([saveName '_max-xAM-Bmode'])
end

function updateROI(src,evt)
    wellIx = str2double(src.Tag(7:end)); % get well index from the tag
    assignin('base', 'wellIx', wellIx);
    
    % pull in variables we need to recompute values
    Nf = evalin('base','Nf');
    Xi = evalin('base','Xi_cell{wellIx}');
    ZixTemp = evalin('base','ZixTemp_cell{wellIx}');
    
    disp(['ROI ' src.Tag ' moved. New position: ' mat2str(evt.CurrentPosition)])
    ROImask = createMask(src); % create mask for ROI
    [z_coords, x_coords] = find(ROImask == 1);
    
    if strcmp(src.Tag(1:5), "sampl")
        % recompute samp ROI and its dimensions
        ZixROI = [min(z_coords) max(z_coords)]; %range of Z coords
        XixROI = [min(x_coords) max(x_coords)]; %range of X coords
        assignin('base', 'ZixROI', ZixROI); %need to assign in base so they can be used by the evalin commands below
        assignin('base', 'XixROI', XixROI);
    elseif strcmp(src.Tag(1:5), "noise")
        % recompute noise ROI and its dimensions
        ZixNoise = [min(z_coords) max(z_coords)]; %range of Z coords
        XixNoise = [min(x_coords) max(x_coords)]; %range of X coords
        assignin('base', 'ZixNoise', ZixNoise); %need to assign in base so they can be used by the evalin commands below
        assignin('base', 'XixNoise', XixNoise);
    end
    
    %update samp rectangle for both the xAM and Bmode figures
    open_figs = findobj('type', 'figure'); %find all figure objects
    for fig_ind = 1:length(open_figs)
        fig = open_figs(fig_ind);
        disp(fig.Tag);
        if strcmp(fig.Tag, 'ROI-layout') %specifically find the layout figure objects
            disp('Found fig');
            roi_axes_array = fig.Children.Children;
            roi_axes = roi_axes_array(length(roi_axes_array) + 1 - wellIx); %get the axes for the updated ROI, need to search in reverse order
            graphics = roi_axes.Children; %get graphics array for specific figure in layout
            for graphic_index = 1:length(graphics) %iterate through graphics array to find samp rectangle
                object = graphics(graphic_index);
                if strcmp(object.Tag, sprintf("sampl_%d", wellIx)) && strcmp(src.Tag(1:5), "sampl")
                    object.Position = [Xi(XixROI(1)) ZixTemp(ZixROI(1)) Xi(XixROI(2))-Xi(XixROI(1)) ZixTemp(ZixROI(2))-ZixTemp(ZixROI(1))]; %update sample ROI position
                elseif strcmp(object.Tag, sprintf("noise_%d", wellIx)) && strcmp(src.Tag(1:5), "noise")
                    object.Position = [Xi(XixNoise(1)) ZixTemp(ZixNoise(1)) Xi(XixNoise(2))-Xi(XixNoise(1)) ZixTemp(ZixNoise(2))-ZixTemp(ZixNoise(1))]; % update noise ROI position
                end
            end
        end
    end
    
    % recalulate sample mean or noise mean and noise STD using the new ROI for all frames
    for frame = 1:Nf
        for imMode = 1:2
            %need to use evalin to run the following commands in the main workspace
            evalin('base',sprintf('ImTemp = Imi{%d,%d,wellIx}(Ixz_cell{wellIx},:);',frame,imMode)); % get each image
            if strcmp(src.Tag(1:5), "sampl")
                evalin('base',sprintf('quants.sampROI_means(%d,%d,wellIx) = mean(mean(ImTemp(ZixROI(1):ZixROI(2),XixROI(1):XixROI(2))));', frame,imMode));
            elseif strcmp(src.Tag(1:5), "noise")
                evalin('base',sprintf('quants.noiseROI_means(%d,%d,wellIx) = mean(mean(ImTemp(ZixNoise(1):ZixNoise(2),XixNoise(1):XixNoise(2))));', frame,imMode));
                evalin('base',sprintf('quants.noiseROI_stds(%d,%d,wellIx) = std2(ImTemp(ZixNoise(1):ZixNoise(2),XixNoise(1):XixNoise(2)));', frame,imMode));
            end
        end
    end
    evalin('base','quantify(quants, VmaxIX);'); % update quants
    
    disp(['Updated ' src.Tag ' mean']);
    
    %update ZixROI and XixROI arrays
    evalin('base','ZixROI_array(wellIx, :) = ZixROI;');
    evalin('base','XixROI_array(wellIx, :) = XixROI;');
    
    % update plots
    evalin('base','make_plots(quants, PlateSize, Nf, confScore,saveName)');
end

function layout = layoutfigures(figs_array,n_rows,n_cols, title_str, cmap, color_lims)
%This function creates the interactive layout of images and the rectangle
%ROIs. The ROIs can be changed, and the callbacks for thise figures will
%still be active

    % figs_array is a gobjects array for the figures that you'd like
    % included in the montage
    f = msgbox('Plotting ROIs...');
    h = figure;
    h.WindowState = 'maximized';
    title(title_str);
    set(h, 'tag', 'ROI-layout');
    layout = tiledlayout(n_rows,n_cols);
    for wellIx = 1:length(figs_array)
        nexttile % create empty axes object and place it in next empty tile of tiled layout
        curr_axes = gca;
        fig = figs_array(wellIx).CurrentAxes; % pull fig object from array
        graphics = fig.Children; %CurrentObject
        for graph_index = 1:length(graphics) %copy over all of the objects in the figure
            old_object = graphics(graph_index);
            new_object = copyobj(old_object,curr_axes,'legacy'); %copy the object into current axes with legacy settings
            if strcmp(old_object.Tag, sprintf("sampl_%d", wellIx)) || strcmp(old_object.Tag, sprintf("noise_%d", wellIx)) %objects with these tags are ROIs and need callbacks
                addlistener(new_object,'ROIMoved',@updateROI); %updateROI is the callback function
            else
                colormap(cmap) % set colormap for images because this isn't copied
                caxis(color_lims) % change color scaling for images because this isn't copied
            end
        end
        row = floor((wellIx-1)/12) + 65;
        col = mod(wellIx-1,12) + 1;
        well = [char(row), num2str(col,'%02d')];
        title(well, 'FontSize', 8);
        set(curr_axes, 'YDir','reverse') % flip image upside-down
        set(curr_axes, 'tag', ['ROI_' num2str(wellIx)]);
        axis off;
    end
    layout.TileSpacing = 'tight';
    layout.Padding = 'tight';
    close(f)
end