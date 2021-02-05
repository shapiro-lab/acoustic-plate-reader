org = params.Stages.Origin;

try;
    sub_Stage_Move_To(params, org);
catch 
    params_old = params;
    
    sub_Close_All_Connections;

    params = sub_AllSettings('BackToOrigin');
    params = sub_Stage_Initialize(params);
    sub_Stage_Move_To(params, org);
end



