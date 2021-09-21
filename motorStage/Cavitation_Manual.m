sub_Close_All_Connections;
params = sub_AllSettings('Cavitation');
params.Scope.averaging = 1;
params.Scope.channels = 4;

params = sub_Scope_Initialize(params);
params = sub_SG_Initialize(params);
params = sub_Stage_Initialize(params);

fprintf(params.Scope.visaObj, ':TRIG:SOURCE EXT');
fprintf(params.Scope.visaObj, ':TRIG:HOLD 9')    
fprintf(params.Scope.visaObj, ':MEAS:CLE');

t_PD = 0.1;

params.SG.Waveform.ch = 1;

params.SG.Waveform.frequency = 4.94E+05;
params.SG.Waveform.voltage = 49.97 * 10^(-params.Amplifier.GainDB/20);

params.SG.Waveform.cycles = params.SG.Waveform.frequency * t_PD;
params.SG.Waveform.period = 1;

close all;
figure('Color', 'w', 'Units', 'inches', 'Position', [0 0 4 6]);

params.Scope.Settings.Position = t_PD / 2;
params.Scope.Settings.TimeRange = t_PD;
params = sub_Scope_ApplySettings(params);
    
params.SG.Waveform.cycles = params.SG.Waveform.frequency * (t_PD);
params.SG.Waveform.period = 1;

disp('MANUAL MODE - saving SG settings only');
% Don't change any SG settings, assume user inputted correctly (finger's
% crossed!)
params.SG.Ch1Status = query(params.SG.visaObj, 'C1:BTWV?');
params.SG.Ch2Status = query(params.SG.visaObj, 'C2:BTWV?');

pause(1);
disp('High fidelity oscilloscope read');
params.Scope.ArmTrigger = 1;
params = sub_Scope_Readout_HQ(params); % Because ArmTrigger = 1, this will cause the Tabor to trigger
disp('Analyzing Data');
params.Results.waveform.XData = params.Scope.Ch(4).XData;
params.Results.waveform.YData = params.Scope.Ch(4).YData;

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

params.Results.w = w(i_plot);
params.Results.t_sample = t_sample;
params.Results.FFT = zeros(n, numel(i_plot));

for i = 0:n-1
    indices = (1 + i * tn_sub) : ((i+1)* tn_sub);
    
    A_sub = A(indices);
    Aw = abs(fft(A_sub));

    params.Results.FFT(i+1,:) = Aw(i_plot);
        
end

 subplot(311);
 plot(t * 1000,A);
 xlabel('Time (ms)')
 ylabel('PCD signal');
 ylim([-0.03 0.03]);

 subplot(312);
 imagesc(params.Results.t_sample * 1000,params.Results.w/1e6,20*log10(params.Results.FFT')); 
 ylabel('Freq (MHz)'); 
 xlabel('Time (ms)');
 c = colorbar;
 c.Label.String = 'dB';
 c.Location = 'NorthOutside';
 set(gca,'CLim',[-50 10])
     
 subplot(313); plot(params.Results.t_sample * 1000, sum(20*log10(params.Results.FFT),2));
 xlabel('Time (ms)');
 ylabel('Sum of dB');
     
 drawnow();

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
sub_Close_All_Connections;

%% Save Results
fld = ['results\' datestr(params.Time, 'yyyy-mm-dd')];
mkdir(fld);

ttl = [params.Name '_' datestr(params.Time, 'yyyy-mm-dd_HH-MM-SS')];

fname = [fld '\results_' ttl '.mat'];
disp(ttl);

save(fname, 'params', '-v7.3')
iname = [fld '\results_' ttl '.png'];
subplot(311); 
title(ttl, 'Interpreter','none');
print(1, '-dpng', iname);
disp('Done!');

clearvars;
