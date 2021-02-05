clc;
sub_Close_All_Connections;
params = sub_AllSettings('MoveSweep');

params = sub_Stage_Initialize(params);

params.Stages.Speed = 100;

loc0 = params.Stages.Origin;

loc1 = loc0;
loc1(params.Stages.y_motor) = loc1(params.Stages.y_motor) - 0.006 / params.Stages.step_distance;

i = 0;

h_tic_master = tic;

while toc(h_tic_master) < 60*5;
    
    if i == 0
        params = sub_Stage_Move_To(params,loc0);
        i = 1;
    else
        params = sub_Stage_Move_To(params,loc1);
        i = 0;
    end
    disp(sprintf('Time elapsed: %1.2f sec', toc(h_tic_master)))
end

params = sub_Stage_Move_To(params,loc0);
disp('Done')