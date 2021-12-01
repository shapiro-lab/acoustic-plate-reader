savedata = 1;
savefig = 0;

% ROI dimensions: voltage, mode, wells
sampROIx = squeeze(sampROI(:,1,:)); % mean intensity of xAM, linear scale
sampROIb = squeeze(sampROI(:,2,:)); % mean intensity of Bmode, linear scale
sampCNRx = squeeze(sampCNR(:,1,:)); % CNR of xAM, dB scale
sampCNRb = squeeze(sampCNR(:,2,:)); % CNR of Bmode, dB scale
noiseROIx = squeeze(noiseROI_mean(:,1,:)); % mean intensity of xAM noise, linear scale
noiseROIb = squeeze(noiseROI_mean(:,2,:)); % mean intensity of Bmode noise, linear scale

PreV = 1:(length(voltage)-1); %Define pre-collapse voltage range
PostV = PreV(end)+1;%Define post-collapse voltage range
split = 3; %number of ramps
PreV_split = reshape(PreV,3,[]);
% sampCNRx(sampCNRx < -7) = -7;

% normx = 20*log10(abs((sampROIx(PreV,:) - noiseROIx(PreV,:)) ./...
%         (sampROIb(PreV,:) - sampROIb(PostV,:)))); % xAM/Bmode, dB scale

% normx = 20*log10(abs((sampROIx(PreV,:) - noiseROIx(PreV,:)) ./...
%         (sampROIb(1,:) - sampROIb(PostV,:)))); % xAM/Bmode, dB scale, normalize only to the first Bmode frame

normx = (sampROIx(PreV,:)); % only plotting raw CNR

vx = PreV;%voltage(PreV); % voltages for plotting

threshold_x = [1:total_n]; % groups/thresholding 
NpS = 4;
numcondition = 5;
groups = cell(1,numcondition);
group_names = strings(1,numcondition);
groups{1} = (1:NpS);
colors_group = linspecer(numcondition);
for gi = 2:numcondition
    groups{gi} = groups{gi-1}(end)+1:groups{gi-1}(end)+NpS;
end
for i = 1:numcondition
    group_names(i) = num2str(i);
end

mnormx_group = [];
normxstd_group = [];

for gi = 1:numcondition
    mnormx_group = [mnormx_group mean(normx(:,groups{gi}),2)];
    groups_data{gi} = normx(:,groups{gi});
    normxstd_group = [normxstd_group std(normx(:,groups{gi}),0,2)/sqrt(length(groups{gi}))];     
end

PlateCoordinate_selected = PlateCoordinate(threshold_x);
selected = normx(:,threshold_x);

colors = linspecer(length(PlateCoordinate_selected));

figure; hold on;
% for i = 1:length(PlateCoordinate_selected)
%     plot(vx,selected(:,i),'LineWidth',2,'DisplayName',PlateCoordinate_selected(i),'Color',colors(i,:));
% end
for i = 1:numcondition
    errorbar(vx,(mnormx_group(1:length(vx),i)),(normxstd_group(1:length(vx),i)),'LineWidth',3,'Color',colors_group(i,:),'DisplayName',group_names(i));            
end
xlabel('Voltage (V)');
legend;
set(gca,'fontsize',16);

sampROIx = squeeze(sampROI(:,1,:)); % mean intensity of xAM, linear scale
sampROIb = squeeze(sampROI(:,2,:)); % mean intensity of Bmode, linear scale
sampCNRx = squeeze(sampCNR(:,1,:)); % CNR of xAM, dB scale
sampCNRb = squeeze(sampCNR(:,2,:)); % CNR of Bmode, dB scale
noiseROIx = squeeze(noiseROI_mean(:,1,:)); % mean intensity of xAM noise, linear scale
noiseROIb = squeeze(noiseROI_mean(:,2,:)); % mean intensity of Bmode noise, linear scale

if savedata
    save([saveName 'data_' datestr(now,'yymmdd-hh-MM-ss')],'sampROI','sampCNR','noiseROI_mean','noiseROI_std','voltage','groups', ...
            'normx','mnormx_group','normxstd_group');
end
if savefig
    savefig([saveName 'normx_' datestr(now,'yymmdd-hh-MM-ss') '.fig'])
end
