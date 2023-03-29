function [hill_funciton_genta,hill_funciton_azithro]=hillEffect(concGen,concAz,parasGenta,parasAzithro)

phi_max=0.67;
phi_min=parasGenta(1);
MIC=parasGenta(3);
Hill_k=parasGenta(2);
hill_funciton_genta=((phi_max-phi_min).*(concGen/MIC).^Hill_k./((concGen/MIC).^Hill_k-(phi_min/phi_max)));


phi_min=parasAzithro(1);
MIC=parasAzithro(3);
Hill_k=parasAzithro(2);
hill_funciton_azithro=((phi_max-phi_min).*(concAz/MIC).^Hill_k./((concAz/MIC).^Hill_k-(phi_min/phi_max)));
end