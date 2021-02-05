% Motor Stage - must be set near the focal point
%
% Function Generator must be manually set at center frequency
% - Run on repeat with timer at 1 ms, for 10 cycles

sub_Close_All_Connections;
params = sub_AllSettings('MoveFromB2toA1');

%% Setup the Initial Motor Stages Parameters
params = sub_Stage_Initialize(params);

%% Search

sub_Stage_Move(params, params.Stages.x_motor, -params.Plate.welldistance / params.Stages.step_distance);
sub_Stage_Move(params, params.Stages.z_motor, -params.Plate.welldistance / params.Stages.step_distance);
        