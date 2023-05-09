sub_Close_All_Connections;
params = sub_AllSettings('Cavitation_Seg');
params.Scope.averaging = 1;
params.Scope.channels = 4;

params = sub_Scope_Initialize(params);
params = sub_SG_Initialize(params);
params = sub_Stage_Initialize(params);

fprintf(params.Scope.visaObj, ':TRIG:SOURCE EXT');
fprintf(params.Scope.visaObj, ':TRIGger:HOLDoff 4e-8')    
fprintf(params.Scope.visaObj, ':MEAS:CLE');

t_PD = 0.01;
segs = 10;

params.SG.Waveform.ch = 1;

params.SG.Waveform.frequency = 4.94E+05;
params.SG.Waveform.voltage = 49.97 * 10^(-params.Amplifier.GainDB/20);

params.SG.Waveform.cycles = params.SG.Waveform.frequency * t_PD;
params.SG.Waveform.period = t_PD * 10;
params.SG.Waveform.repeats = segs;

close all;
figure('Color', 'w', 'Units', 'inches', 'Position', [0 0 4 6]); clf;

params.Scope.Settings.Position = t_PD / 2;
params.Scope.Settings.TimeRange = t_PD;
params = sub_Scope_ApplySettings(params);

fprintf(params.Scope.visaObj, ':ACQuire:MODE SEGMented')
fprintf(params.Scope.visaObj, ':ACQuire:SEGMented:COUNt %1.0f', segs)

disp('Configuring SG to send single pulse');
params = sub_SG_Stop(params);
params = sub_SG_ApplySettingsForTrigger(params);
params = sub_SG_Start(params);
pause(1);
disp('High fidelity oscilloscope read');
params.Scope.ArmTrigger = 1;
params = sub_Scope_Readout_HQ(params); % Because ArmTrigger = 1, this will cause the Tabor to trigger
disp('Analyzing Data');

params.Results.waveform.XData = params.Scope.Ch(4).XData;
params.Results.waveform.YData = params.Scope.Ch(4).YData;

dt_sample = 0.001;

t = params.Results.waveform.XData;
A = params.Results.waveform.YData; 
A = A - mean(A(:));

subplot(311);
plot(t * 1000,A, 'k-');
xlabel('Time (ms)')
ylabel('PCD signal');
ylim([-0.03 0.03]);

for seg = 1:segs
t = params.Results.waveform.XData(:,seg);
A = params.Results.waveform.YData(:,seg); 
A = A - mean(A);
    
t_sample = min(t):dt_sample:max(t)-dt_sample;
n = numel(t_sample);

tn_sub = floor(numel(t) / n);
fs = 1/(t(2)-t(1)); %Sampling frequency
fft_pts = tn_sub; % Nb points
w = (0:fft_pts-1)./fft_pts.*fs;
i_plot = 50:floor(numel(w)/2);

w = w(i_plot);
FFT = [];

for i = 0:n-1
    indices = (1 + i * tn_sub) : ((i+1)* tn_sub);
    
    A_sub = A(indices);
    Aw = abs(fft(A_sub));

    FFT(i+1,:) = Aw(i_plot);
        
end



 subplot(312);
 imagesc(t_sample * 1000,w/1e6,20*log10(FFT')); hold on;
 ylabel('Freq (MHz)'); 
 xlabel('Time (ms)');
 xlim([0, 1000*t_PD * 10 *segs]);
 set(gca,'CLim',[-50 10])

     
 subplot(313); plot(t_sample * 1000, sum(20*log10(FFT),2), 'k-'); hold on;
 xlabel('Time (ms)');
 ylabel('Sum of dB');
     
 drawnow();
end

subplot(312);
 c = colorbar;
 c.Label.String = 'dB';
 c.Location = 'NorthOutside';
 


disp('Saving data');

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
params = sub_SG_Stop(params);
sub_Close_All_Connections;

%% Save Results
fld = ['results\' datestr(params.Time, 'yyyy-mm-dd')];
mkdir(fld);

ttl = [params.Name '_' datestr(params.Time, 'yyyy-mm-dd_HH-MM-SS')];

fname = [fld '\results_' ttl '.mat'];
disp(ttl);
disp(sprintf('Energy = %1.2e', sum(params.Results.waveform.YData(:).^2)));
disp(sprintf('Energy = %+1.2f dB', 10*log10(sum(params.Results.waveform.YData(:).^2))));

params = rmfield(params, 'Scope');

save(fname, 'params', '-v7.3')
iname = [fld '\results_' ttl '.png'];
subplot(311); 
title(ttl, 'Interpreter','none');
print(1, '-dpng', iname);
disp('Done!');

clearvars;
