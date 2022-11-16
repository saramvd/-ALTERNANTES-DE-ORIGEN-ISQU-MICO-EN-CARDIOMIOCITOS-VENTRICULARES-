clc
%SOFTWARE PARA LA REALIZACIÓN DE SIMULACIONES ESTÁTICAS DE FORMA SISTEMÁTICA

%Inicializamos las variables del modelo
corr_act.IKATP=0;
corr_act.INaK=0;
corr_act.SERCA=0;
corr_act.IpCa=0;
corr_act.INa=0;
corr_act.INaL=0;
corr_act.ICaL=0;
corr_act.NCX=0;

batchYN = 1;

isq_act.O2=1;
isq_act.pH=1;
isq_act.Met=1;
isq_act.wash_out_isq = 1;

settings.minIsqTotal = 10.0;
settings.minIsqIni = 0.0;

%Variaciones en los valores de los metabolitos
%Primera simulación
% ATP_var = [7.4 7.733 8.066 8.4];
% ADP_var = [316 320 324 328];
% pHi_var = [6.65 6.7166 6.7833 6.85];
% LPC_var = [4.5 5 5.5 6];
% Ko_var = [10.5 10.933 11.366 11.8];

%Segunda simulación
ATP_var = [6 7.33 8.66 10];
ADP_var = [100 183.33 266.66 350];
pHi_var = [6.2 6.533 6.866 7.2];
LPC_var = [2 4 6 8];
Ko_var = [11.5 11.833 12.166 12.5];

%Bucle de las simulaciones
for idx_ATP=1:length(ATP_var)
    for idx_ADP=1:length(ADP_var)
        for idx_pHi=1:length(pHi_var)
            for idx_LPC=1:length(LPC_var)
                for idx_Ko=1:length(Ko_var)
                    
                    %definimos los valores que tendrán los settings en cada bucle
                    settings.ATPi_ini = ATP_var(idx_ATP);
                    settings.ADPi_ini = ADP_var(idx_ADP);
                    settings.pHi_ini = pHi_var(idx_pHi);
                    settings.pHo_ini = pHi_var(idx_pHi)+0.2;
                    settings.LPC_ini = LPC_var(idx_LPC);
                    settings.ko_ini = Ko_var(idx_Ko);                    
                   
                    %Llamamos al modelo
                    [minTi,Ti,StateVars,currents]=main2019_ORd_MMChAL(settings,isq_act,corr_act,batchYN);
                    
                    %Ponemos nombre a cada una de las simulaciones
                    simulacion_minTi = strcat('minTi_ATP',string(ATP_var(idx_ATP)),'_ADP',string(ADP_var(idx_ADP)),'_pHi',string(pHi_var(idx_pHi)),'_LPC',string(LPC_var(idx_LPC)),'_ko',string(Ko_var(idx_Ko)),'.mat');
                    simulacion_StateVars = strcat('StateVars_ATP',string(ATP_var(idx_ATP)),'_ADP',string(ADP_var(idx_ADP)),'_pHi',string(pHi_var(idx_pHi)),'_LPC',string(LPC_var(idx_LPC)),'_ko',string(Ko_var(idx_Ko)),'.mat');
                    simulacion_currents = strcat('currents_ATP',string(ATP_var(idx_ATP)),'_ADP',string(ADP_var(idx_ADP)),'_pHi',string(pHi_var(idx_pHi)),'_LPC',string(LPC_var(idx_LPC)),'_ko',string(Ko_var(idx_Ko)),'.mat');
                    
                    %Guardamos los vectores de cada simulación
                    save(simulacion_minTi, 'minTi');
                    save(simulacion_StateVars, 'StateVars');
                    save(simulacion_currents, 'currents');
                  
                    %Figura del Potencial de acción (los últimos cinco)
                    plot(minTi,StateVars(:,1));
                    xlim([4.9, 5]); %establecemos valor del eje x para que solo sean los ultimos cinco PA
                    simulacion_PA = strcat('PA ATP-',string(ATP_var(idx_ATP)),' ADP-',string(ADP_var(idx_ADP)),' pHi-',string(pHi_var(idx_pHi)),' LPC-',string(LPC_var(idx_LPC)),' ko-',string(Ko_var(idx_Ko)));
                    title(simulacion_PA);
                    
                    nombre_figura = strcat('FiguraPA_ATP',string(ATP_var(idx_ATP)),'_ADP',string(ADP_var(idx_ADP)),'_pHi',string(pHi_var(idx_pHi)),'_LPC',string(LPC_var(idx_LPC)),'_ko',string(Ko_var(idx_Ko)),'.fig');
                    saveas (gcf, nombre_figura);
                    
                end
            end
        end
    end
end 
