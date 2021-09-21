params = sub_AllSettings('SetHSTestingParameters');

sub_Close_All_Connections;
params = sub_SG_Initialize(params);

disp(query(params.SG.visaObj, 'C1: OUTP?'))
disp(query(params.SG.visaObj, 'C1: BSWV?'))
disp(query(params.SG.visaObj, 'C1: BTWV?'))
disp(query(params.SG.visaObj, 'C2: OUTP?'))
disp(query(params.SG.visaObj, 'C2: BSWV?'))
disp(query(params.SG.visaObj, 'C2: BTWV?'))