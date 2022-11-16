function [minTi,Ti,StateVars,currents]=main2019_ORd_MMChAL(settings,isq_act,corr_act,batchYN)    

%[minTi,Ti,StateVars,currents]=main_ORd_MMChA(settings,isq_act,corr_act,batchYN);
%===================================================================
% MATLAB Implementation of the O'Hara-Rudy dynamic (ORd) model
% modified by Maite, Patri, Ana & Chema
%===================================================================
% SETTINGS
% settings.minIsqIni: La isquemia empieza al inicio de este minuto (defecto: 5)
% settings.minIsqTotal: La isquemia dura estos minutos (defecto: 30)
% settings.BCL: pacing cycle length in ms (defecto: 1000)
% settings.ATPi_ini: [ATP]i inicial en mM (defecto: 10)
% settings.ATPi_fin: [ATP]i final en mM (defecto: 2.5)
% settings.ADP_type: Tipo de variación gráfica del ADP (defecto: 3) 1:Weiss; 2:Clarke;3: Smith
% settings.ADPi_ini: [ADP]i inicial en microM (defecto: 15)
% settings.ADPi_fin_1: [ADP]i final en microM tras la primera fase de ascenso (defecto: 80)
% settings.ADPi_fin_2: [ADP]i final en microM tras la segunda fase de ascenso o descenso (defecto: 80)
% settings.pHi_ini: pHi inicial (defecto: 7.2)
% settings.pHi_fin: pHi final (defecto: 6.2)
% settings.pHo_ini: pHo inicial (defecto: 7.4)
% settings.pHo_fin: pHo final (defecto: 6.4)
% settings.LPC_ini: valor inicial de LPC en microM (defecto: 2)
% settings.LPC_fin: valor final de LPC en microM (defecto: 15)
% settings.ko_ini: valor inicial de Ko en mM (defecto: 5.4)
% settings.minDT_LPC: minutos durante los cuales evoluciona la LPC (defecto: 30)
% settings.wash_out_norm: 0 si no hay wash-out en normoxia, 1 (defecto) si hay wash-out en normoxia
% settings.tau_wash_out: tau de wash out en ms(defecto: 15000)
%===================================================================
% ISQ_ACT: '1' para activar el componente isqémico; '0' para desactivarlo
% isq_act.O2: afecta al ATPi y al ADPi
% isq_act.pH: afecta a la acidosis intracelular y extracelular
% isq_act.Met: afecta a los metabolitos (p.ej. LPC)
% isq_act.wash_out_isq: 0 (defecto) si no hay wash-out en isquemia, 1 si hay wash-out en isquemia
%===================================================================
% CORR_ACT: '1' para activar el efecto isquémico sobre una corriente
% corr_act.IKATP
% corr_act.INaK
% corr_act.SERCA
% corr_act.IpCa
% corr_act.INa
% corr_act.INaL
% corr_act.ICaL
% corr_act.NCX
%===================================================================
%batchYN: '1' is simulation is in batch mode (many sequentially programmed)

%===================================================================
% minTi: time vector in minutos
% Ti: time vector in ms
% StateVars: state vars matrix
% currents: structure with currents and other stuff
%===================================================================



settings = setDefaultSettings(settings);    %Parámetros de configuración
isq_act = setDefaultIsqAct(isq_act);    % Activación/desactivación de componentes isquémicos
corr_act = setDefaultCorrAct(corr_act); % Activación/desactivación de efectos isquémicos sobre corrientes


CL = settings.BCL;
minSimTotal = settings.minIsqIni + settings.minIsqTotal; % min totales de la simulación 
msSimTotal = minSimTotal * 60 * 1000; % ms totales de la simulación
msIsqIni = settings.minIsqIni * 60 * 1000;   % La isquemia empieza en este milisegundo tras el inicio de la simulación
msIsqTotal = settings.minIsqTotal * 60 * 1000;   % La isquemia empieza en este milisegundo tras el inicio de la simulación
beats = round(msSimTotal/CL)-1;  %beats=1080%Lucia +240;

%initial conditions for state variables. There are 42 state variables. 
%Para una única célula en un estado inicial de 1Hz(diástole)
v=-87; %potencial de reposo mV
nai=7; %Concentración Na intracelular (todas las [] están en mM.
nass=nai;
ki=145;
kss=ki;
cai=1.0e-4;
cass=cai;
cansr=1.2;
cajsr=cansr;
m=0;
hf=1;
hs=1;
j=1;
hsp=1;
jp=1;
mL=0;
hL=1;
hLp=1;
a=0;
iF=1;
iS=1;
ap=0;
iFp=1;
iSp=1;
d=0;
ff=1;
fs=1;
fcaf=1;
fcas=1;
jca=1;
nca=0;
ffp=1;
fcafp=1;
xrf=0;
xrs=0;
xs1=0;
xs2=0;
xk1=1;
Jrelnp=0;
Jrelp=0;
CaMKt=0;

%X0 is the vector for initial sconditions for state variables
X0=[v nai nass ki kss cai cass cansr cajsr m hf hs j hsp jp mL hL hLp a iF iS ap iFp iSp d ff fs fcaf fcas jca nca ffp fcafp xrf xrs xs1 xs2 xk1 Jrelnp Jrelp CaMKt]';

options=[];%options for ode solver
%options= odeset('RelTol',1e-7);
StateVars = [];
Ti=[];
Currents = [];

%=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0
t_start = tic;
progress=0;
t0 = clock;
if batchYN==0
    h4 = waitbar(0,'1','Name',' Calculando Potencial de Acción Nº           ','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    setappdata(h4,'canceling',0);
end
%=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0=0

for n=[1:beats];
    [time X]=ode15s(@model2019_ORd_MMChAL,[0 CL],X0,options,settings,isq_act,corr_act,msIsqIni,msIsqTotal,n,1);

    if batchYN==0
         % Check for clicked Cancel button
         if getappdata(h4,'canceling')
            break
        end
   % Update waitbar and message
        waitbar(n/beats,h4,sprintf('%12.0f',n))
    end
    
    X0=X(size(X,1),:);    
    Ti=[Ti; time+(CL*(n-1))];
    StateVars=[StateVars; X];   
end
if batchYN==0
    delete(h4)
end




minTi_proto = Ti/(60*1000);   %vector de tiempo en minutos
minTi = -settings.minIsqIni + ((settings.minIsqTotal + settings.minIsqIni) / (max(minTi_proto))) * (minTi_proto);
%rename values in the state variables vector
v=X(:,1);
nai=X(:,2);
nass=X(:,3);
ki=X(:,4);
kss=X(:,5);
cai=X(:,6);
cass=X(:,7);
cansr=X(:,8);
cajsr=X(:,9);
m=X(:,10);
hf=X(:,11);
hs=X(:,12);
j=X(:,13);
hsp=X(:,14);
jp=X(:,15);
mL=X(:,16);
hL=X(:,17);
hLp=X(:,18);
a=X(:,19);
iF=X(:,20);
iS=X(:,21);
ap=X(:,22);
iFp=X(:,23);
iSp=X(:,24);
d=X(:,25);
ff=X(:,26);
fs=X(:,27);
fcaf=X(:,28);
fcas=X(:,29);
jca=X(:,30);
nca=X(:,31);
ffp=X(:,32);
fcafp=X(:,33);
xrf=X(:,34);
xrs=X(:,35);
xs1=X(:,36);
xs2=X(:,37);
xk1=X(:,38);
Jrelnp=X(:,39);
Jrelp=X(:,40);
CaMKt=X(:,41);


%calculate and name dependent variables (i.e. currents and fluxes)

if batchYN==0
    h5 = waitbar(1,' Computando corrientes... ');
end

for i=1:size(Ti);
    beat_actual=floor(Ti(i)/CL)+1;
    IsJs=model2019_ORd_MMChAL(Ti(i),StateVars(i,:),settings,isq_act,corr_act,msIsqIni,msIsqTotal,beat_actual,0);
    currents.INa(i)=IsJs(1);
    currents.INaL(i)=IsJs(2);
    currents.Ito(i)=IsJs(3);
    currents.ICaL(i)=IsJs(4);
    currents.IKr(i)=IsJs(5);
    currents.IKs(i)=IsJs(6);
    currents.IK1(i)=IsJs(7);
    currents.INaCa_i(i)=IsJs(8);
    currents.INaCa_ss(i)=IsJs(9);
    currents.INaK(i)=IsJs(10);
    currents.IKb(i)=IsJs(11);
    currents.INab(i)=IsJs(12);
    currents.ICab(i)=IsJs(13);
    currents.IpCa(i)=IsJs(14);
    currents.Jdiff(i)=IsJs(15);
    currents.JdiffNa(i)=IsJs(16);
    currents.JdiffK(i)=IsJs(17);
    currents.Jup(i)=IsJs(18);
    currents.Jleak(i)=IsJs(19);
    currents.Jtr(i)=IsJs(20);
    currents.Jrel(i)=IsJs(21);
    currents.CaMKa(i)=IsJs(22);
    currents.Istim(i)=IsJs(23);
    currents.IKATP(i)=IsJs(24);
    currents.ICaK(i)=IsJs(25);
    currents.f_ATP_KATP(i) = IsJs(26);
    currents.ATPi(i) = IsJs(27);
    currents.ADPi(i) = IsJs(28);
    currents.EK(i) = IsJs(29);
    currents.pHi(i) = IsJs(30);
    currents.pHo(i) = IsJs(31);
    currents.Km(i) = IsJs(32);
    currents.HATP(i) = IsJs(33);
    currents.f_pHo_INa(i)=IsJs(34);
    currents.f_LPC_INaL(i)=IsJs(35);
    currents.f_pHi_ICaL(i)=IsJs(36);
    currents.f_pHo_ICaL(i)=IsJs(37);
    currents.f_pHi_NCX(i)=IsJs(38);
    currents.LPC(i)=IsJs(39);
    currents.f_LPC_INa(i)=IsJs(40);
    currents.f_ATP_INaK(i)=IsJs(41);
    currents.f_ATP_IpCa(i)=IsJs(42);
    currents.f_ATP_Iup(i)=IsJs(43); 
    currents.f_pHi_INaK(i)=IsJs(44); 
    currents.wash_out(i)=IsJs(45); 
    currents.fINap(i)=IsJs(46); 
    
    currents.f_pHi_Jup(i)=IsJs(47);
    currents.f_pHo_NCX(i)=IsJs(48);
    currents.f_pHo_INaL(i)=IsJs(49);
    currents.f_IKr_pHo(i)=IsJs(50);
    currents.f_IKr_LPC(i)=IsJs(51);
    currents.f_pHi_INa(i)=IsJs(52);
    
    currents.ko(i)=IsJs(53);
    
%     currents.INa(i)=IsJs(1);
%     currents.INaL(i)=IsJs(2);
%     currents.Ito(i)=IsJs(3);
%     currents.ICaL(i)=IsJs(4);
%     currents.IKr(i)=IsJs(5);
%     currents.IKs(i)=IsJs(6);
%     currents.IK1(i)=IsJs(7);
%     currents.INaCa_i(i)=IsJs(8);
%     currents.INaCa_ss(i)=IsJs(9);
%     currents.INaK(i)=IsJs(10);
%     currents.IKb(i)=IsJs(11);
%     currents.INab(i)=IsJs(12);
%     currents.ICab(i)=IsJs(13);
%     currents.IpCa(i)=IsJs(14);
%     currents.Jdiff(i)=IsJs(15);
%     currents.JdiffNa(i)=IsJs(16);
%     currents.JdiffK(i)=IsJs(17);
%     currents.Jup(i)=IsJs(18);
%     currents.Jleak(i)=IsJs(19);
%     currents.Jtr(i)=IsJs(20);
%     currents.Jrel(i)=IsJs(21);
%     currents.CaMKa(i)=IsJs(22);
%     currents.Istim(i)=IsJs(23);
%     currents.IKATP(i)=IsJs(24);
%     currents.ICaK(i)=IsJs(25);
%     currents.f_ATP_KATP(i) = IsJs(26);
%     currents.ATPi(i) = IsJs(27);
%     currents.ADPi(i) = IsJs(28);
%     currents.EK(i) = IsJs(29);
%     currents.pHi(i) = IsJs(30);
%     currents.pHo(i) = IsJs(31);
%     currents.Km(i) = IsJs(32);
%     currents.HATP(i) = IsJs(33);
%     currents.f_pHi_INa(i)=IsJs(34);
%     currents.f_LPC_INaL(i)=IsJs(35);
%     currents.f_pHi_ICaL(i)=IsJs(36);
%     currents.f_pHo_ICaL(i)=IsJs(37);
%     currents.f_pHi_NCX(i)=IsJs(38);
%     currents.LPC(i)=IsJs(39);
%     currents.f_LPC_INa(i)=IsJs(40);
%     currents.f_ATP_INaK(i)=IsJs(41);
%     currents.f_ATP_IpCa(i)=IsJs(42);
%     currents.f_ATP_Iup(i)=IsJs(43); 
%     currents.f_pHi_INaK(i)=IsJs(44); 
%     currents.wash_out(i)=IsJs(45); 
    
    %[time IsJs]=model_ORdmm,[0 CL],IsJs,CL,n,options,0);
    %ISJS=[IsJs(end,1:end)]; 
    %Currents=[Currents;IsJs];
    %waitbar(i/size(Ti),h5)
   % waitbar(beat_actual/beats,h5)

end
if batchYN==0
    delete(h5)
end

if batchYN==0
    t_elapsed = toc(t_start);
else
    t_elapsed = toc(t_start);
end

if batchYN==0
    mensaje = ['Simulation of ' int2str(minTi_proto(end)) ' minutes of biological time lasted for ' int2str(round(t_elapsed)) ' seconds of real time'];
	msgbox(mensaje,'Simulation complete')
end

save('monitor00.mat');

end

function settings = setDefaultSettings(settings)
    if ~isfield(settings, 'minIsqTotal'), settings.minIsqTotal = 5.0; end
    if ~isfield(settings, 'minIsqIni'), settings.minIsqIni = 0.0; end
    if ~isfield(settings, 'BCL'), settings.BCL = 1000; end
    if ~isfield(settings, 'ATPi_ini'), settings.ATPi_ini = 10; end
    if ~isfield(settings, 'ATPi_fin'), settings.ATPi_fin = 1; end
    if ~isfield(settings, 'ADP_type'), settings.ADP_type = 3; end    
    if ~isfield(settings, 'ADPi_ini'), settings.ADPi_ini = 15; end
    if ~isfield(settings, 'ADPi_fin_1'), settings.ADPi_fin_1 = 80; end
    if ~isfield(settings, 'ADPi_fin_2'), settings.ADPi_fin_2 = 80; end
    if ~isfield(settings, 'pHi_ini'), settings.pHi_ini = 7.2; end
    if ~isfield(settings, 'pHi_fin'), settings.pHi_fin = 5.9; end
    if ~isfield(settings, 'pHo_ini'), settings.pHo_ini = 7.4; end
    if ~isfield(settings, 'pHo_fin'), settings.pHo_fin = 6.1; end
    if ~isfield(settings, 'fmINaL'), settings.fmINaL = 8; end
    if ~isfield(settings, 'LPC_ini'), settings.LPC_ini = 2; end
    if ~isfield(settings, 'LPC_fin'), settings.LPC_fin = 20; end
    if ~isfield(settings, 'ko_ini'), settings.ko_ini = 5.4; end
    if ~isfield(settings, 'minDT_LPC'), settings.minDT_LPC = 30; end
    if ~isfield(settings, 'wash_out_norm'), settings.wash_out_norm = 1.0; end
    if ~isfield(settings, 'tau_wash_out'), settings.tau_wash_out = 15000; end
        
    
end

function isq_act = setDefaultIsqAct(isq_act)
    if ~isfield(isq_act, 'O2'), isq_act.O2 = 1; end
    if ~isfield(isq_act, 'pH'), isq_act.pH = 1; end
    if ~isfield(isq_act, 'Met'), isq_act.Met = 1; end
    if ~isfield(isq_act, 'wash_out_isq'), isq_act.wash_out_isq = 0.0; end    
end

function corr_act = setDefaultCorrAct(corr_act)
    if ~isfield(corr_act, 'IKATP'), corr_act.IKATP = 1; end
    if ~isfield(corr_act, 'INaK'), corr_act.INaK = 1; end
    if ~isfield(corr_act, 'SERCA'), corr_act.SERCA = 1; end
    if ~isfield(corr_act, 'IpCa'), corr_act.IpCa = 1; end
    if ~isfield(corr_act, 'INa'), corr_act.INa = 1; end
    if ~isfield(corr_act, 'INaL'), corr_act.INaL = 1; end
    if ~isfield(corr_act, 'ICaL'), corr_act.ICaL = 1; end
    if ~isfield(corr_act, 'NCX'), corr_act.NCX = 1; end
end




% v: 1
% nai: 2
% nass: 3
% ki: 4
% kss: 5
% cai: 6
% cass: 7
% cansr: 8
% cajsr: 9
% m: 10
% hf: 11
% hs: 12
% j: 13
% hsp: 14
% jp: 15
% mL: 16
% hL: 17
% hLp: 18
% a: 19
% iF: 20
% iS: 21
% ap: 22
% iFp: 23
% iSp: 24
% d: 25
% ff: 26
% fs: 27
% fcaf: 28
% fcas: 29
% jca: 30
% nca: 31
% ffp: 32
% fcafp: 33
% xrf: 34
% xrs: 35
% xs1: 36
% xs2: 37
% xk1: 38
% Jrelnp: 39
% Jrelp: 40
% CaMKt: 41

