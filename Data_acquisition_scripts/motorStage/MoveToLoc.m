function params = MoveToLoc(desiredloc)

try; delete(instrfind); end

params = sub_Stage_Initialize(struct);
params = sub_Stage_Update_Positions(params);
params.Stages.MoveDelay = 1;
params = sub_Stage_Move_To(params, desiredloc);

end