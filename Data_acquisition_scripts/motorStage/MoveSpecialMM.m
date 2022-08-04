sub_Close_All_Connections;

%% Setup the Initial Motor Stages Parameters
params = sub_AllSettings('MoveSpecial');
params = sub_Stage_Initialize(params);

%For Wavelength chamber, each well is separated by 0.5"
% v2 ds_mm = 0.5 * 25.4;
% 37.08 mm in newest model

% For 670 kHz transudcer 0.072 m, 48 us;
% For 500 khz transducer 0.0325, 21.6 us;
% For 100 kHz transducer 0.05268;


%% Search

dx_mm = (34.2 - 32.5);
dy_mm = 0;
dz_mm = 0;

sub_Stage_Move(params, params.Stages.z_motor, (dz_mm / 1000) / params.Stages.step_distance);
sub_Stage_Move(params, params.Stages.y_motor, (dy_mm / 1000) / params.Stages.step_distance);
sub_Stage_Move(params, params.Stages.x_motor, (dx_mm / 1000) / params.Stages.step_distance);