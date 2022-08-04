params = sub_AllSettings('SetTestingParameters');

sub_Close_All_Connections;
params = sub_SG_Initialize(params);


params.SG.Waveform.ch = 1;
params.SG.Waveform.frequency = 4.94E+05;
params.SG.Waveform.voltage = 49.97 * 10^(-params.Amplifier.GainDB/20);
cycles = 0.1 * params.SG.Waveform.frequency;
params.SG.Waveform.period = cycles / params.SG.Waveform.frequency * 10;
params.SG.Waveform.cycles = cycles;

params = sub_SG_ApplySettings(params);
params = sub_SG_Stop(params);

params = sub_Stage_Initialize(params); 

ax = params.Stages.x_motor;
org = params.Stages.Origin;
step = 37.08e-3 / params.Stages.step_distance;
duration = 120;

loc = zeros(3,11);

loc(:,1) = org;
for i = 2:11
    loc(:,i) = org;
    loc(ax,i) = loc(ax,i) + (i-1) * step;
end

for i = 1:4
   disp(sprintf('Moving to Position %1.0f', i));
   params = sub_Stage_Move_To(params, loc(:,i));
   params = sub_SG_Start(params);
   tic;
   
   disp_amp = 8e-3 / params.Stages.step_distance;
   
   while toc < duration
       for dy = -disp_amp:disp_amp/10:disp_amp
           loc2 = loc(:,i);
           loc2(params.Stages.y_motor) = loc2(params.Stages.y_motor) + dy;
           parmas = sub_Stage_Move_To(params, loc2);
           pause(0.1);
       end
   end
   params = sub_SG_Stop(params);
end

params = sub_Stage_Move_To(params, org);