sub_Close_All_Connections;
clearvars;

params = sub_AllSettings('Cavitation_MultiLoc');
params.Scope.averaging = 1;
params.Scope.channels = 4;

params = sub_Scope_Initialize(params);
params = sub_SG_Initialize(params);
params = sub_Stage_Initialize(params);

fld = ['results\' datestr(params.Time, 'yyyy-mm-dd')];
mkdir(fld);
ttl = [params.Name '_' datestr(params.Time, 'yyyy-mm-dd_HH-MM-SS')];
fname = [fld '\results_' ttl '.mat'];
iname = [fname(1:end-4) '.gif'];
i1name = [fname(1:end-4) '.png'];
ilocname = [fname(1:end-4) '_Loc.png'];

fprintf(params.Scope.visaObj, ':TRIG:HOLD 9')
fprintf(params.Scope.visaObj, ':MEAS:CLE');

t_PD = 0.1;

close all;
figure('Color', 'w', 'Units', 'inches', 'Position', [0 0 4 6]);

params.SG.Waveform.ch = 1;
params.SG.Waveform.frequency = 4.94E+05;
params.SG.Waveform.voltage = 49.97 * 10^(-params.Amplifier.GainDB/20);

params.SG.Waveform.cycles = params.SG.Waveform.frequency * t_PD;
params.SG.Waveform.period = t_PD * 10;

params.Scope.Settings.Position = t_PD / 2;
params.Scope.Settings.TimeRange = t_PD;
params = sub_Scope_ApplySettings(params);

wavelength = params.Acoustic.MediumAcousticSpeed / params.SG.Waveform.frequency;
dim1 = params.Stages.z_motor;
step_ds = (wavelength / 24) / params.Stages.step_distance;
step_n = 24;
org = params.Stages.Origin;
locs = zeros(3,step_n);
for i = 1:step_n
    ds = (i - step_n/2) * step_ds;
    locs(:,i) = org;
    locs(dim1,i) = org(dim1) + ds;
end
params.Scan.Locations = locs;

h_tic = tic;
params = sub_SG_Stop(params);
params = sub_SG_ApplySettingsForTrigger(params);
params = sub_SG_Start(params);
pause(3);




for j = 1:120

    r = 1+floor(rand*step_n);
    params.Results.LocationIndex(j) = r;
    params = sub_Stage_Move_To(params, params.Scan.Locations(:,r));

    params.Scope.ArmTrigger = 1;
    params = sub_Scope_Readout_HQ(params); % Because ArmTrigger = 1, this will cause the Tabor to trigger   
    params.Results.TimeOfAq(j) = toc(h_tic);
    
    disp(sprintf('%1.0f - %1.2f s', j, params.Results.TimeOfAq(j)));
    
    params.Results.waveforms{j} = sub_Data_CompressWaveform(params.Scope.Ch(4));
    
	t = params.Scope.Ch(4).XData;
    A = params.Scope.Ch(4).YData; 
    A = A - mean(A);
    
    dt_sample = 0.001;

    t_sample = 0:dt_sample:t_PD-dt_sample;
    n = numel(t_sample);

    tn_sub = floor(numel(t) / n);
    fs = 1/(t(2)-t(1)); %Sampling frequency
    fft_pts = tn_sub; % Nb points
    w = (0:fft_pts-1)./fft_pts.*fs;
    i_plot = 50:floor(numel(w)/2);
    
    w = w(i_plot);
        
    data = zeros(n, numel(i_plot));

    for i = 0:n-1
        indices = (1 + i * tn_sub) : ((i+1)* tn_sub);

        A_sub = A(indices);
        Aw = abs(fft(A_sub));

        data(i+1,:) = Aw(i_plot);

    end
    
    params.Results.WaveEnergy(j) = sum(A.^2);
    
     figure(2); subplot(211);
     hold off;
     for i = 1:j
        s = params.Scan.Locations(dim1,params.Results.LocationIndex(i)) - params.Stages.Origin(dim1);
        s = 1000 * params.Stages.step_distance * s;
       
        scatter(s, 10*log10(params.Results.WaveEnergy(i)), 25, [1-i/120,0,0])
     hold on;
     end
     
     xlabel('Location (mm)');
     ylabel('PCD Signal Energy (dB)');
     
     subplot(212);
     hold off;
     for i = 1:j
     scatter(i, 10*log10(params.Results.WaveEnergy(i)), 25, [1-i/120,0,0])
     hold on;
     end
     
     xlabel('Pulse Number');
     ylabel('PCD Signal Energy (dB)');
     
    
     figure(1);
     subplot(311);
     plot(t * 1000,A);
     xlabel('Time (ms)')
     ylabel('PCD signal');
     ylim([-0.03 0.03]);

     subplot(312);
     imagesc(t_sample * 1000, w/1e6, 20*log10(data')); 
     ylabel('Freq (MHz)'); 
     xlabel('Time (ms)');
     c = colorbar;
     c.Label.String = 'dB';
     c.Location = 'NorthOutside';
     set(gca,'CLim',[-50 10])

     subplot(313); plot(t_sample * 1000, sum(20*log10(data),2));
     xlabel('Time (ms)');
     ylabel('Sum of dB');

     drawnow();
    
           frame = getframe(gcf); 
      im = frame2im(frame); 
      [imind, cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if j == 1 
          imwrite(imind,cm,iname,'gif', 'Loopcount',1,'DelayTime',1/3); 
          imwrite(imind,cm,i1name,'png'); 
      else 
          imwrite(imind,cm,iname,'gif','WriteMode','append','DelayTime',1/3); 
      end 
     
end

figure(2); 
frame = getframe(gcf); 
im = frame2im(frame); 
[imind, cm] = rgb2ind(im,256); 
imwrite(imind,cm,ilocname,'png'); 

%% Save the code form this script in params file
try
    IO = fopen([mfilename '.m'],'r');
    params.ScriptCode = textscan(IO,'%s','Delimiter','\n'); 
    fclose(IO);
    clear IO;
catch
    disp('Could not save a copy of script code in params file')
end


%% Close All Connections
params = sub_Stage_Move_To(params, params.Stages.Origin);
sub_Close_All_Connections;

%% Save Results

params = rmfield(params, 'Scope');

summary = struct();
summary.locs = locs;
summary.LocationIndex = params.Results.LocationIndex;
summary.WaveEnergy = params.Results.WaveEnergy;

a = whos('params');
clear A t data A_sub Aw indices i_plot w t_sample frame im imind cm

disp(sprintf('- Saving Data (%1.0f MB)', a.bytes/(1e6)));
disp(['- ' ttl]);
disp('- Please wait...');
save(fname, '-v7.3')
disp('- Done!');

clearvars
