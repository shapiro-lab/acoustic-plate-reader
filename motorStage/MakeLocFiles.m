files = dir('\results\2018-03-27\*Loc*.mat');

params.Stages.step_distance = 0.0254/10/400; 

close all;
figure('Color', 'w', 'Units', 'inches', 'Position', [0 0 4 6]);

for ii = 1:numel(files);
   load(['\results\2018-03-27\' files(ii).name], 'summary');
   
   dim1 = 0;
   
   for i = 1:3
       if summary.locs(i,1) ~= summary.locs(i,2)
           dim1 = i;
       end
   end
   
    subplot(211);
     hold off;
     for i = 1:numel(summary.WaveEnergy)
        s = summary.locs(dim1,summary.LocationIndex(i)) - median(summary.locs(dim1,:));
        s = 1000 * params.Stages.step_distance * s;
       
        scatter(s, 10*log10(summary.WaveEnergy(i)), 25, [1-i/120,0,0])
     hold on;
     end
     
     ylim([-10 15]);
     xlabel('Location (mm)');
     ylabel('PCD Signal Energy (dB)');
     
     subplot(212);
     hold off;
     for i = 1:numel(summary.WaveEnergy)
     scatter(i, 10*log10(summary.WaveEnergy(i)), 25, [1-i/120,0,0])
     hold on;
     end
     
     xlabel('Pulse Number');
     ylabel('PCD Signal Energy (dB)');
     ylim([-10 15]);
     
     f = files(ii).name;
     f = [f(1:end-4) '_Loc.png'];
     
     
        subplot(211);
     title(f, 'Interpreter','none');
     
      frame = getframe(gcf); 
      im = frame2im(frame); 
      [imind, cm] = rgb2ind(im,256); 
      imwrite(imind,cm,f,'png'); 
end