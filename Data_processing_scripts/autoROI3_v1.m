% clear all
% close all
% 
% pathName = 'G:\.shortcut-targets-by-id\0B24ONICaZ0z9djczVE1ZR3BnWU0\Shapiro Lab Information\Data\Rob\96-well_plate_scans\GvpA-B-mutants\A-lib-2\A-lib-K22R-A2\A-lib-K22R-A2_P3_4_stable-37C_P_R1_C3';

% PlateProc1_v2
% 
close all
clear iZt ixZtemp
%%
sample_depth = [1 10]; % display depth range in mm
ROI_depth = [3 8]; % sample depth range in mm

%Mohamed 12/27: DisplayMode is now obsolete, only thing that shows is layout and 96 well plots
DisplayMode = 1; % How to display images: 0 = show nothing, 1 = always show images + ROIs, 2 = only show images + ROIs with low-confidence ROI selection

testfilter = 0; % 0 or 1 to display the thresholding filter
xROI_size = 10; % size of the ROI in x dimension    
zROI_size = 40; % size of the ROI in z dimension 
filter_size1 = [40 5]; % size of the median filter in [z,x] dimension
noise_stdratio = 3; % ratio of the noise STD added to the noise mean to decide noise floor 
noise_slices = repmat([7 8],total_n,1); % default noise depth (in mm) for noise ROI selection
noise_slices_backup = [8 9;9 10;10 11]; % backup noise depth (in mm) for noise ROI selection in case of abnormal noise level
noiseThresh = 100; % threshold of noise level (mV) for abnormal high noise level / changing noise depth slices
TemplateMode = 3; % ROI template to use: 1 = filled wells, 2/3 = empty wells with thinner (2) or thicker (3) well wall, use 1 as default and 2/3 if well is not full
if ~exist('WellTemplates','var') % load the ROI template if needed, please put the .mat file in accessible folders
    load('WellTemplates.mat');
end
WellTemplate = WellTemplates{TemplateMode}; % Apply ROI template as selected.
CST = ConfidenceScoreThreshold(TemplateMode); % Apply corresponding confidence score threshold for reminder to adjust ROIs (abitrary defined thresholds)

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
sampROI_means = nan(Nf,2,total_n); % initialize array for sample means
noiseROI_means = nan(Nf,2,total_n); % initialize array for noise means
noiseROI_stds = nan(Nf,2,total_n); % initialize array for noise ROI STDs
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
            imagesc(Xi, ZixTemp, 20*log10(abs(ImTemp)), [20 80]); axis image; colormap bone;
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
%             ImTemp_xAM = Imi{Nf-1,1,wellIx}; % Fetch xAM image, currently using second-to-last voltage image
            ImTemp_xAM = Imi{10,1,wellIx};
            ImTemp_xAM = ImTemp_xAM(ixZtemp,:);
            imagesc(Xi, ZixTemp, 20*log10(abs(ImTemp_xAM)), [20 80]); axis image; colormap hot;
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
            imagesc(Xi, ZixTemp, 20*log10(abs(ImTemp)), [20 80]); axis image; colormap bone;
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
            imagesc(Xi, ZixTemp, 20*log10(abs(ImTemp_xAM)), [20 80]); axis image; colormap hot;
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
                ImTemp = Imi{frame,imMode,wellIx}(ixZtemp,:);
                sampROI_means(frame,imMode,wellIx) = mean(mean(ImTemp(ZixROI(1):ZixROI(2),XixROI(1):XixROI(2))));
                noiseROI_means(frame,imMode,wellIx) = mean(mean(ImTemp(iZd_noise,:)));
                noiseROI_stds(frame,imMode,wellIx) = std2(ImTemp(iZd_noise,:));
            end
        end
    end
    waitbar(wellIx/total_n)
end
close(h);

sampSBR = sampROI_means ./ noiseROI_means; % Calculate sample SBR
sampSBR_dB = 20 * log10(sampSBR); % Calculate sample SBR in dB
sampSBR_diff = sampSBR(1:10,:,:) - sampSBR(11:20,:,:); % Calculate pre-/post-collapse difference SBR
sampSBR_diff_dB = 20 * log10(abs(sampSBR_diff)); % Calculate pre-/post-collapse difference SBR in dB

sampCNR = (sampROI_means - noiseROI_means) ./ noiseROI_stds; % Calculate sample CNR
sampCNR_dB = 20 * log10(abs(sampCNR)); % Calculate sample CNR in dB
sampCNR_diff = sampCNR(1:10,:,:) - sampCNR(11:20,:,:); % Calculate pre-/post-collapse difference CNR
sampCNR_diff_dB = 20 * log10(abs(sampCNR_diff)); % Calculate pre-/post-collapse difference CBR in dB

AM_Bmode_ratio_dB = 20*log10(abs((sampROI_means(:,1,:) - noiseROI_means(:,1,:)) ./ (sampROI_means(:,2,:) - noiseROI_means(:,1,:)))); % xAM/Bmode, dB scale

% save data
clear pressure imMode;
save([saveName '_data_' datestr(now,'yymmdd-hh-MM-ss')],'P','saveName','sampROI_means','sampCNR','AM_Bmode_ratio_dB','noiseROI_means','noiseROI_stds','voltage','PlateCoordinate','PlateSize','ROI_Centers','confScore','Nf');

%% plot ROI quants with microplateplot
% reshape ROI CNRs
% sampCNRs_dB new dimensions: well rows, well columns, frames, imaging modes
sampCNRs_dB = permute(reshape(sampCNR_dB, Nf, 2, PlateSize(2), PlateSize(1)), [4 3 1 2]);
sampCNRs_diff_dB = permute(reshape(sampCNR_diff_dB, Nf/2, 2, PlateSize(2), PlateSize(1)), [4 3 1 2]);
AM_Bmode_ratios_dB = permute(reshape(AM_Bmode_ratio_dB, Nf, 1, PlateSize(2), PlateSize(1)), [4 3 1 2]);
confScore = permute(reshape(confScore, PlateSize(2), PlateSize(1)), [1 2]);

% find max signal achieved by each sample at any voltage
maxs_AM = squeeze(max(sampCNRs_dB, [], 3));
maxs_AM_diff = squeeze(max(sampCNRs_diff_dB, [], 3));
maxs_AM_Bmode_ratio = squeeze(max(AM_Bmode_ratios_dB, [], 3));

% make and save microplate plots
plot_xAM(maxs_AM, saveName)
plot_AM_Bmode_ratio(maxs_AM_Bmode_ratio, saveName)

Bmode_layout = layoutfigures(Bmode_figs,PlateSize(1),PlateSize(2), 'Bmode layout', 'bone'); %create Bmode figure layout
xAM_layout = layoutfigures(xAM_figs,PlateSize(1),PlateSize(2), 'xAM layout', 'hot'); %create xAM figure layout
%% function definitions
function make_plots(sampCNR, PlateSize, Nf, AM_Bmode_ratio,confScore,saveName)
    % This function makes the maxs_AM and maxs_AM_Bmode_ratio. If they
    % exist, it deletes the previous maxs figures and recomputes the values
    
    %reshape ROI CNRs
    % sampCNR new dimensions: well rows, well columns, frames, imaging modes
    sampCNRs = permute(reshape(sampCNR, Nf, 2, PlateSize(2), PlateSize(1)), [4 3 1 2]);
    AM_Bmode_ratios = permute(reshape(AM_Bmode_ratio, Nf, 1, PlateSize(2), PlateSize(1)), [4 3 1 2]);
    confScore = permute(reshape(confScore, PlateSize(2), PlateSize(1)), [1 2]);
    
    % find max signal achieved by each sample at any voltage
    maxs_AM = squeeze(max(sampCNRs, [], 3));
    maxs_AM_Bmode_ratio = squeeze(max(AM_Bmode_ratios, [], 3));
    
    %delete previous plots
    open_figs = findobj('type', 'figure');
    for fig_ind = 1:length(open_figs)
        fig = open_figs(fig_ind);
        if strcmp(fig.CurrentAxes.Tag, 'max-xAM') || strcmp(fig.CurrentAxes.Tag, 'max-xAM:Bmode')%this will delete the previous maxs figures if they exist, we remake them later
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
    sampMask = createMask(src); % create mask for samp ROI
    
    % recompute samp ROI and its dimensions
    [z_coords, x_coords] = find(sampMask == 1);
    ZixROI = [min(z_coords) max(z_coords)]; %range of Z coords
    XixROI = [min(x_coords) max(x_coords)]; %range of X coords
    assignin('base', 'ZixROI', ZixROI);%need to assign in base so they can be used by the evalin commands below
    assignin('base', 'XixROI', XixROI);
    
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
            for graph_index = 1:length(graphics) %iterate through graphics array to find samp rectangle
                object = graphics(graph_index);
                if strcmp(object.Tag, sprintf("sampl_%d", wellIx))
                    object.Position = [Xi(XixROI(1)) ZixTemp(ZixROI(1)) Xi(XixROI(2))-Xi(XixROI(1)) ZixTemp(ZixROI(2))-ZixTemp(ZixROI(1))];%update the position of the rectangle
                end
            end
        end
    end
    
    % recalulate sample mean, noise mean, and noise STD using the new ROI for all frames
    for frame = 1:Nf
        for imMode = 1:2
            %need to use evalin to run the following commented out commands in the main workspace
            %ImTemp = Imi{frame,imMode,wellIx}(ixZtemp,:);
            % sampROI_means(frame,imMode,wellIx) = mean(mean(ImTemp(ZixROI(1):ZixROI(2),XixROI(1):XixROI(2))));
            %noiseROI_means(frame,imMode,wellIx) = mean(mean(ImTemp(iZd_noise,:)));
            %noiseROI_stds(frame,imMode,wellIx) = std2(ImTemp(iZd_noise,:));
            evalin('base',sprintf('ImTemp = Imi{%d,%d,wellIx}(Ixz_cell{wellIx},:);',frame,imMode));
            evalin('base',sprintf('sampROI_means(%d,%d,wellIx) = mean(mean(ImTemp(ZixROI(1):ZixROI(2),XixROI(1):XixROI(2))));', frame,imMode));
            evalin('base',sprintf('noiseROI_means(%d,%d,wellIx) = mean(mean(ImTemp(iZd_noise,:)));', frame,imMode));
            evalin('base',sprintf('noiseROI_stds(%d,%d,wellIx) = std2(ImTemp(iZd_noise,:));', frame,imMode));
        end
    end
    
    %need to use evalin to run the following commented out commands in the main workspace
    %sampCNR_dB(:,:,wellIx) = 20 * log10(abs(sampROI_means(:,:,wellIx) - noiseROI_means(:,:,wellIx)) ./ noiseROI_stds(:,:,wellIx)); % Replace sample CNR
    %AM_Bmode_ratio_dB = 20*log10(abs((sampROI_means(:,1,:) - noiseROI_means(:,1,:)) ./ (sampROI_means(:,2,:) - noiseROI_means(:,1,:)))); % Recalculate xAM/Bmode, dB scale
    evalin('base','sampCNR_dB(:,:,wellIx) = 20 * log10(abs(sampROI_means(:,:,wellIx) - noiseROI_means(:,:,wellIx)) ./ noiseROI_stds(:,:,wellIx));'); % Replace sample CNR
    evalin('base','AM_Bmode_ratio_dB = 20*log10(abs((sampROI_means(:,1,:) - noiseROI_means(:,1,:)) ./ (sampROI_means(:,2,:) - noiseROI_means(:,1,:))));'); % Recalculate xAM/Bmode, dB scale
    disp(['Updated ' src.Tag ' mean']);
    
    %update ZixROI and XixROI arrays
    evalin('base','ZixROI_array(wellIx, :) = ZixROI;');
    evalin('base','XixROI_array(wellIx, :) = XixROI;');
    
    % update plots
    evalin('base','make_plots(sampCNR_dB, PlateSize, Nf, AM_Bmode_ratio_dB,confScore,saveName)');
 
    % save updated quants
    evalin('base', "save([saveName '_data_' datestr(now,'yymmdd-hh-MM-ss')],'P','saveName','sampROI_means','sampCNR_dB','AM_Bmode_ratio_dB','noiseROI_means','noiseROI_stds','voltage','PlateCoordinate','PlateSize','ROI_Centers','confScore','Nf');")
end

function layout = layoutfigures(figs_array,n_rows,n_cols, title_str, cmap)
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
    colormap(cmap)
    layout = tiledlayout(n_rows,n_cols);
    for wellIx = 1:length(figs_array)
        nexttile
        curr_axes = gca;
        fig = figs_array(wellIx).CurrentAxes; % pull fig object from array
        graphics = fig.Children; %CurrentObject
        for graph_index = 1:length(graphics) %copy over all of the objects in the figure
            old_object = graphics(graph_index);
            new_object = copyobj(old_object,curr_axes,'legacy'); %copy over the object with legacy
            if strcmp(old_object.Tag, sprintf("sampl_%d", wellIx)) || strcmp(old_object.Tag, sprintf("noise_%d", wellIx)) %objects with these tags need callbacks
                addlistener(new_object,'ROIMoved',@updateROI); %updateROI is the callback function
            end
        end
        set(curr_axes, 'YDir','reverse')
        set(curr_axes, 'tag', ['ROI_' num2str(wellIx)]);
        axis off;
    end
    layout.TileSpacing = 'compact';
    layout.Padding = 'compact';
    close(f)
end