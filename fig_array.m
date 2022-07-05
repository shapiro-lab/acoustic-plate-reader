
close all
clear all

total_n = 12;
data = randi(256,1000,1000,2,total_n);

sampROI_means = nan(1,2,total_n);
noiseROI_means = nan(1,2,total_n);
noiseROI_stds = nan(1,2,total_n);

figs = gobjects(1,2,total_n); % make array to hold figure objects

for wellIx = 1:total_n

        figure(wellIx);
        img = image(data(:,:,2,wellIx)); % use Bmode image 

        hold on;
        noiseROI = drawrectangle('Position',[100,100,100,100],'EdgeColor','w','LineWidth',2, 'Tag',['noise_' num2str(wellIx)]);
        addlistener(noiseROI,'ROIMoved',@updateNoiseROI);
        sampROI = drawrectangle('Position',[200,200,100,100],'EdgeColor','g','LineWidth',2, 'Tag',['samp_' num2str(wellIx)]);
        addlistener(sampROI,'ROIMoved',@updateSampROI);
        title(['Frame #' num2str(wellIx)]);
        figs(wellIx) = figure(wellIx); % store figure object in array

        % Create mask for ROIs
        noiseMask = createMask(noiseROI);
        sampMask = createMask(sampROI);

        % Multiply images by masks and take means of masked ROIs
        noiseROI_means(wellIx) = mean(nonzeros(data.*noiseMask), 'all', 'omitnan');
        sampROI_means(wellIx) = mean(nonzeros(data.*sampMask), 'all', 'omitnan');
end


% sampROI_means(frame,imMode,wellIx) = mean(mean(ImTemp(ZixROI(1):ZixROI(2),XixROI(1):XixROI(2))));
% noiseROI_means(frame,imMode,wellIx) = mean(mean(ImTemp(iZd_noise,:)));
% noiseROI_stds(frame,imMode,wellIx) = std2(ImTemp(iZd_noise,:));


function updateNoiseROI(src,evt)
    disp(['ROI ' src.Tag ' moved. New position: ' mat2str(evt.CurrentPosition)])
    noiseMask = createMask(src); % Create mask for noise ROI
    
    % Calculate new quants
    wellIx = str2double(src.Tag(7:end));
    meanTemp = mean(nonzeros(evalin('base', sprintf('data(:,:,:, %d)', wellIx)).*noiseMask), [1 2 4], 'omitnan')
    stdTemp = std2(nonzeros(evalin('base', sprintf('data(:,:,:, %d)', wellIx)).*noiseMask));
    
    % Assign new quants into array
    evalin('base',sprintf('noiseROI_means(1,:, %d)=%d;',wellIx,meanTemp));
    disp(['New ' src.Tag ' mean: ' num2str(meanTemp)]);
    evalin('base',sprintf('noiseROI_stds(1,:, %d)=%d;',wellIx,stdTemp));
    disp(['New ' src.Tag ' STD: ' num2str(stdTemp)]);
end

function updateSampROI(src,evt)
    disp(['ROI ' src.Tag ' moved. New position: ' mat2str(evt.CurrentPosition)])
    sampMask = createMask(src); % Create mask for samp ROI
    
    % Calculate new quant
    wellIx = str2double(src.Tag(6:end));
    meanTemp = mean(nonzeros(evalin('base', sprintf('data(:,:, %d)', wellIx)).*sampMask), 'all', 'omitnan');
    
    % Assign new quant into array
    evalin('base',sprintf('sampROI_means(1, %d)=%d;',wellIx,meanTemp));
    disp(['New ' src.Tag ' mean: ' num2str(meanTemp)]);
end