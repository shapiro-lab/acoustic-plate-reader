
close all
num_wells = 12;
figs = gobjects(1,num_wells);
noiseQuants = zeros(1,num_wells);
sampQuants = zeros(1,num_wells);
data = randi(256,1000,1000,num_wells);

for figIdx = 1:num_wells
        
        figure(figIdx);
        img = image(data(:,:,figIdx));
        
        hold on;
        noiseROI = drawrectangle('Position',[100,100,100,100],'EdgeColor','w','LineWidth',2, 'Tag',['noise_' num2str(figIdx)]);
        addlistener(noiseROI,'ROIMoved',@updateNoiseROI);
        sampROI = drawrectangle('Position',[200,200,100,100],'EdgeColor','g','LineWidth',2, 'Tag',['samp_' num2str(figIdx)]);
        addlistener(sampROI,'ROIMoved',@updateSampROI);
        title(['Frame #' num2str(figIdx)]);
        figs(figIdx) = figure(figIdx); %store figure object in array
        
        % Create mask for ROIs
        noiseMask = createMask(noiseROI);
        sampMask = createMask(sampROI);
        
        % Multiply images by masks and take means of masked ROIs
        noiseQuants(figIdx) = mean(nonzeros(data.*noiseMask), 'all', 'omitnan');
        sampQuants(figIdx) = mean(nonzeros(data.*sampMask), 'all', 'omitnan');
end





function updateNoiseROI(src,evt)
    disp(['ROI ' src.Tag ' moved. New position: ' mat2str(evt.CurrentPosition)])
    noiseMask = createMask(src); % Create mask for noise ROI
    
    % Calculate new quant
    figIdx = str2double(src.Tag(7:end));
    temp = mean(nonzeros(evalin('base', sprintf('data(:,:, %d)', figIdx)).*noiseMask), 'all', 'omitnan');
    
    % Assign new quant into array
    evalin('base',sprintf('noiseQuants(1, %d)=%d;',figIdx,temp));
    disp(['New quant: ' num2str(temp)]);
end

function updateSampROI(src,evt)
    disp(['ROI ' src.Tag ' moved. New position: ' mat2str(evt.CurrentPosition)])
    sampMask = createMask(src); % Create mask for samp ROI
    
    % Calculate new quant
    figIdx = str2double(src.Tag(6:end));
    temp = mean(nonzeros(evalin('base', sprintf('data(:,:, %d)', figIdx)).*sampMask), 'all', 'omitnan');
    
    % Assign new quant into array
    evalin('base',sprintf('sampQuants(1, %d)=%d;',figIdx,temp));
    disp(['New quant: ' num2str(temp)]);
end