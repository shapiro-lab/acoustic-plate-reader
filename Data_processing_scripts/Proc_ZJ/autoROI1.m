%overlay = 1;
showf5 = 0;

%Ind = 59;
Imt = Im{1,2,1};
%dImi(:,:,imgMode,Ind);
filter_size1 = [20 5];
filter_size2 = [5 5];
filter_sigma = 2.5;
filter_box = [9 1];
%SizeThreshold = 300;

threshold_ratio = 1;
auto_noise = 1;
noiseSize = [2 2];
sample_depth = [2 7];
sample_width = 10;
iZt = find(Zi>=sample_depth(1)):find(Zi>=sample_depth(2));
filter_depth = zeros(size(Imt));
filter_depth(iZt,:) = 1;

if auto_noise
    dZ = mean(diff(Zi));
    dX = mean(diff(Xi));
    noiseTest = [(length(Zi) - round(noiseSize(1)/dZ)) (length(Xi) - round(noiseSize(2)/dX)) ];
    noiseROIt = [mean(mean(Imt(noiseTest(1):end,noiseTest(2):end))) std2(Imt(noiseTest(1):end,noiseTest(2):end))];
else
    noiseROIt = noiseROI(:,Ind);
end

Imt_f1 = medfilt2(Imt,filter_size1);
%threshold1 = (threshold_ratio * (noiseROIt(1) + 2 * noiseROIt(2)));
if imgMode == 1
    threshold1 = 10;
elseif imgMode == 2
    threshold1 = 0.1;
end

Imt_f2 = (Imt_f1 > threshold1);
Imt_f3 = Imt_f2 .* filter_depth;
%Imt_f41 = medfilt2(Imt_f3,filter_size2);
Imt_f4 = imgaussfilt(Imt_f3,filter_sigma);
% Imt_f42 = double(Imt_f41>=1);
% Imt_f4 = imboxfilt(Imt_f42,filter_box);
Imt_f5 = (Imt_f4>=1);
if overlay
    fig = figure('Position',[877 1 561 800]);
    ax1 = axes;
    hb = imagesc(Xi, Zi, 20*log10(abs(Imt)), [20 80]);axis image;colormap hot
    if overlay == 1
        ax2 =axes;
        hx = imagesc(Xi, Zi,Imt_f5,[0 1]); axis image
        linkaxes([ax1,ax2])
        ax2.Visible = 'off';
        ax2.XTick = [];
        ax2.YTick = [];
        colormap(ax2,'cool')
        alphadata = Imt_f5* 0.5;
        set(hx,'AlphaData',alphadata);
    end
elseif showf5
    fig = figure('Position',[877 1 561 800]);
    imagesc(Xi, Zi,Imt_f5,[0 1]); axis image
    colormap hot
end
ROISize = sum(sum(Imt_f5));
ROISize(ROISize<1) = 1;
autoROI = sum(sum(Imt .* Imt_f5));

% stats = regionprops('table',Imt_f5,'Centroid','BoundingBox');
% if height(stats) ~= 1
%     %sandp = [sandp i];
%     disp('Warning: Number of centroid ~= 1. Auto detection may fail.')
% else
% %     imagesc(Xi,Zi,Imt,[20 80]), axis image
% %     colormap hot
% %     ROI_center = round(stats.Centroid);
% %     sample_width = min(sample_width,stats.BoundingBox(3)*dX);
% %     rectangle('Position',[Xi(ROI_center(1)) - sample_width/2 Zi(iZt(1)) sample_width (Zi(iZt(end)) - Zi(iZt(1)))],...
% %         'EdgeColor','w','LineWidth',2);
% end