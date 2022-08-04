%%
close all

savedata = 1;
savefigure = 1;
groupalong = 'columns';
split = 1; % number of voltage ramps
NpS = 1; % number of replicates
numcondition = 96; % number of unique samples

%%
% ROI dimensions: voltage, mode, wells
sampROIx = squeeze(sampROI(:,1,:)); % mean intensity of xAM, linear scale
sampROIb = squeeze(sampROI(:,2,:)); % mean intensity of Bmode, linear scale
sampCNRx = squeeze(sampCNR(:,1,:)); % CNR of xAM, dB scale
sampCNRb = squeeze(sampCNR(:,2,:)); % CNR of Bmode, dB scale
noiseROIx = squeeze(noiseROI_mean(:,1,:)); % mean intensity of xAM noise, linear scale
noiseROIb = squeeze(noiseROI_mean(:,2,:)); % mean intensity of Bmode noise, linear scale

PreV = 1:42; % Define pre-collapse voltage range
% PostV = PreV(end)+1; % Define post-collapse voltage
PostV = length(voltage)-PreV(end)+1:length(voltage); % Define post-collapse voltage range
PreV_split = reshape(PreV,[],split)'; % split pre-collapse voltage range into individual ramps

% normx = 20*log10(abs((sampROIx(PreV,:) - noiseROIx(PreV,:)) ./...
%         (sampROIb(PreV,:) - sampROIb(PostV,:)))); % xAM/Bmode, dB scale
% yLabel = 'xAM/Bmode CNR (dB)';

% normx = 20*log10(abs((sampROIx(PreV,:) - noiseROIx(PreV,:)) ./...
%         (sampROIb(1,:) - sampROIb(PostV,:)))); % xAM/Bmode, dB scale, normalize only to the first Bmode frame
% yLabel = 'xAM/Bmode CNR (dB), normalized only to the first Bmode frame';

normx = (sampCNRx(:,:)); % only plot raw CNR
yLabel = 'xAM CNR (dB)';

% normx = 20*log10(abs((sampROIx(PreV,:) - sampROIx(PostV,:)))); % xAM/Bmode, dB scale
% yLabel = 'xAM pre- post-collapse (dB)';

vx = voltage(PreV_split(1,:)); %voltage(PreV); % voltages for plotting
%% group samples
% threshold_x = [1:total_n]; % groups/thresholding 
groups = cell(1,numcondition);
group_names = strings(1,numcondition);
colors_group = linspecer(numcondition);

if strcmp(groupalong, 'rows')
    groups{1} = (1:NpS);
    for gi = 2:numcondition
        groups{gi} = groups{gi-1}(end) + 1:groups{gi-1}(end) + NpS;
    end
    for sample = 1:numcondition
        group_names(sample) = num2str(sample);
    end
else
    num_stripe = [(numcondition-mod(numcondition,P.zLines))/P.zLines,mod(numcondition,P.zLines)];
    for gj = 1:num_stripe(1)
        for gi = 1:P.zLines
            groups{gi + P.zLines * (gj - 1)} = (gi + NpS*(gj-1)*P.zLines):P.zLines:((NpS-1)*P.zLines + gi + NpS*(gj-1)*P.zLines);
            group_names(gi + P.zLines * (gj - 1)) = num2str(gi + P.zLines * (gj - 1));
        end
    end
    if num_stripe(2)
        for gi = 1:num_stripe(2)
            groups{gi + P.zLines * num_stripe(1)} = (gi + NpS*num_stripe(1)*P.zLines):P.zLines:((NpS-1)*P.zLines + gi + NpS*num_stripe(1)*P.zLines);
            group_names(gi + P.zLines * num_stripe(1)) = num2str(gi + P.zLines * num_stripe(1));
        end
    end
end

mnormx_group = [];
normxstd_group = [];

for gi = 1:numcondition
    mnormx_group = [mnormx_group mean(normx(:,groups{gi}),2)];
    groups_data{gi} = normx(:,groups{gi});
    normxstd_group = [normxstd_group std(normx(:,groups{gi}),0,2)/sqrt(length(groups{gi}))];     
end

% PlateCoordinate_selected = PlateCoordinate(threshold_x);
% selected = normx(:,threshold_x);

%% plot voltage ramps
figure; hold on;
% for i = 1:length(PlateCoordinate_selected)
%     plot(vx,selected(:,i),'LineWidth',2,'DisplayName',PlateCoordinate_selected(i),'Color',colors(i,:));
% end

% [1:15; 16:32; 33:51; 52:72]
% ind = 5:15;
% ind = 20:32;
% ind = 37:51;
% ind = 56:72;
% vx = vx(ind);
% PreV_split = PreV_split(ind);
for samp = 1:numcondition
    for ramp = 1:split
        plot(vx,mnormx_group(PreV_split(ramp,:),samp),'LineWidth',3,'Color',[colors_group(samp,:) 1 - 0.3*(ramp-1)],'DisplayName',group_names(samp));
%         pl = errorbar(vx,(mnormx_group(PreV_split(j,:),i)),(normxstd_group(PreV_split(j,:),i)),'LineWidth',3,'Color',[colors_group(i,:) 1 - 0.3*(j-1)],'DisplayName',group_names(i));
    end
end
xlabel('Voltage (V)');
ylabel(yLabel)
legend;
set(gca,'fontsize',16);

if savefigure
    savefig([saveName '_normx_' datestr(now,'yymmdd-hh-MM-ss') '.fig'])
end

%%
sampROIx = squeeze(sampROI(:,1,:)); % mean intensity of xAM, linear scale
sampROIb = squeeze(sampROI(:,2,:)); % mean intensity of Bmode, linear scale
sampCNRx = squeeze(sampCNR(:,1,:)); % CNR of xAM, dB scale
sampCNRb = squeeze(sampCNR(:,2,:)); % CNR of Bmode, dB scale
noiseROIx = squeeze(noiseROI_mean(:,1,:)); % mean intensity of xAM noise, linear scale
noiseROIb = squeeze(noiseROI_mean(:,2,:)); % mean intensity of Bmode noise, linear scale

% if savedata
%     save([saveName 'data_' datestr(now,'yymmdd-hh-MM-ss')],'sampROI','sampCNR','noiseROI_mean','noiseROI_std','voltage','groups', ...
%             'normx','mnormx_group','normxstd_group');
% end
