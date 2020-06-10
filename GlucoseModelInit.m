%clear all
%clc
%warning off

%% Havorka Average Glucose-Insulin Model Constants
%{
k_12 = .066; %(1/min)  Transfer rate [Havorka et al 2002]
k_a1 = .006; %(1/min) Deactivation rate [Havorka et al 2002]
k_a2 = .06; %(1/min) Deactivation rate [Havorka et al 2002]
k_a3 = .03; %(1/min) Deactivation rate [Havorka et al 2002]
k_e = .138; %(1/min) Insulin elimination from plasma [Havorka et al 1993]
V_I = .12; %(L/kg) Insulin distribution volume [Havorka et al 1993]
V_G = .16; %(L/kg) Glucose distribution volume [Havorka et al 2002]
Bio = .8; %(unitless) Carbohydrate Bioavailability [Livesey et al 1998]
t_max = 40; %(min) Time-to-maximum of CHO absorption [Livesey et al 1998]
U_G_ceil = .0275 %(mmol/kg/min) Maximum gut glucose flux drawn from U [.02, .035]
R_cl = .003;
R_thr = 9;

% Calculated Average Model Parameters
S_IT = 1/(51.2 * 10^-4); %(L/min/mU) Insulin sensitivity of distribution/transport [Havorka et al 2002]
S_ID = 1/(8.2 *10^-4); %(L/min/mU) Insulin sensitivity of disposal [Havorka et al 2002]
S_IE = 1/(520 * 10^-4); %(L/min/mU) Insulin sensitivity of EGP [Havorka et al 2002]
EGP_0 = .0161; %(mmol/kg/min) EGP extrapolated to zero insulin concentration [Havorka et al 2002]
Fs_01 = .0097/8.5; %(mmol/kg/min) Non-insulin-dependent glucose flux [Havorka et al 2002]
t_maxI = 55; %(min) Time-to-maximum of absorption of subcutaneously injected short-acting insulin [Howey et al 1994], [Rave et al 1999]

k_b1 = k_a1 / S_IT;
k_b2 = k_a2 / S_ID;
k_b3 = k_a3 / S_IE;
%}

%% Havorka Patient 1 Specific Parameters

kg = 69; %(kg) Patient Weight
needs = 14; %(U/day) patient insulin needs

Fs_01 = 12.95 * 10^-3 * kg; %(mmol/kg/min) Non-insulin-dependent glucose flux [Havorka et al 2002]
EGP_0 = 19.62 * 10^-3 * kg; %(mmol/kg/min) EGP extrapolated to zero insulin concentration [Havorka et al 2002]
k_12 = 10.95 * 10^-2; %(1/min)  Transfer rate [Havorka et al 2002]
V_G = 17.87 * 10^-2 * kg; %(L/kg) Glucose distribution volume [Havorka et al 2002]
R_thr = 11.70; %(mmol/L)
R_cl = 1.19 * 10^-2; %(1/min)

S_IT = 77.1 * 10^-4; %(L/min/mU) Insulin sensitivity of distribution/transport [Havorka et al 2002]
S_ID = 3.14 * 10^-4; %(L/min/mU) Insulin sensitivity of disposal [Havorka et al 2002]
S_IE = 377 * 10^-4; %(L/min/mU) Insulin sensitivity of EGP [Havorka et al 2002]
k_b1 = .21 * 10^-2;
k_b2 = 39.56 * 10^-2;
k_b3 = 8.03 * 10^-2;
k_a1 = S_IT * k_b1; %(1/min) Deactivation rate [Havorka et al 2002]
k_a2 = S_ID * k_b2; %(1/min) Deactivation rate [Havorka et al 2002]
k_a3 = S_IE * k_b3; %(1/min) Deactivation rate [Havorka et al 2002]

k_a = 1.98 * 10^-2; %(1/min)
k_e = 13.2 * 10^-2; %(1/min) Insulin elimination from plasma [Havorka et al 1993]
V_I = 11.3 * 10^-2 * kg; %(L/kg) Insulin distribution volume [Havorka et al 1993]

Bio = .71; %(unitless) Carbohydrate Bioavailability [Livesey et al 1998]
t_max = 43; %(min) Time-to-maximum of CHO absorption [Livesey et al 1998]
U_G_ceil = .0275 * kg; %(mmol/kg/min) Maximum gut glucose flux drawn from U [.02, .035]

k_a_int = 7.36 * 10^-2;

t_maxI = 55; %(min) Time-to-maximum of absorption of subcutaneously injected short-acting insulin [Howey et al 1994], [Rave et al 1999]


G_0 = 15.3;

%% Havorka Patient 2 Specific Parameters
%{
kg = 90; %(kg) Patient Weight
needs = 39; %(U/day) patient insulin needs

Fs_01 = 9.77 * 10^-3 * kg; %(mmol/kg/min) Non-insulin-dependent glucose flux [Havorka et al 2002]
EGP_0 = 9.10 * 10^-3 * kg; %(mmol/kg/min) EGP extrapolated to zero insulin concentration [Havorka et al 2002]
k_12 = 5.09 * 10^-2; %(1/min)  Transfer rate [Havorka et al 2002]
V_G = 14.50 * 10^-2 * kg; %(L/kg) Glucose distribution volume [Havorka et al 2002]
R_thr = 9.22; %(mmol/L)
R_cl = 1.30 * 10^-2; %(1/min)

S_IT = 11.0 * 10^-4; %(L/min/mU) Insulin sensitivity of distribution/transport [Havorka et al 2002]
S_ID = 1.58 * 10^-4; %(L/min/mU) Insulin sensitivity of disposal [Havorka et al 2002]
S_IE = 73 * 10^-4; %(L/min/mU) Insulin sensitivity of EGP [Havorka et al 2002]
k_b1 = .06 * 10^-2;
k_b2 = 1.36 * 10^-2;
k_b3 = 2.02 * 10^-2;
k_a1 = S_IT * k_b1; %(1/min) Deactivation rate [Havorka et al 2002]
k_a2 = S_ID * k_b2; %(1/min) Deactivation rate [Havorka et al 2002]
k_a3 = S_IE * k_b3; %(1/min) Deactivation rate [Havorka et al 2002]

out = [k_b1*10^2 k_b2*10^2 k_b3*10^2 k_a1 k_a2 k_a3]

k_a = 1.60 * 10^-2; %(1/min)
k_e = 10.1 * 10^-2; %(1/min) Insulin elimination from plasma [Havorka et al 1993]
V_I = 13.1 * 10^-2 * kg; %(L/kg) Insulin distribution volume [Havorka et al 1993]

Bio = .90; %(unitless) Carbohydrate Bioavailability [Livesey et al 1998]
t_max = 55; %(min) Time-to-maximum of CHO absorption [Livesey et al 1998]
U_G_ceil = .0275 * kg; %(mmol/kg/min) Maximum gut glucose flux drawn from U [.02, .035]

k_a_int = 15.10 * 10^-2;

t_maxI = 55; %(min) Time-to-maximum of absorption of subcutaneously injected short-acting insulin [Howey et al 1994], [Rave et al 1999]
%}

%% Havorka Patient 4 Specific Parameters
%{
kg = 87; %(kg) Patient Weight
needs = 23; %(U/day) patient insulin needs

Fs_01 = 8.09 * 10^-3 * kg; %(mmol/kg/min) Non-insulin-dependent glucose flux [Havorka et al 2002]
EGP_0 = 8.25 * 10^-3 * kg; %(mmol/kg/min) EGP extrapolated to zero insulin concentration [Havorka et al 2002]
k_12 = 6.35 * 10^-2; %(1/min)  Transfer rate [Havorka et al 2002]
V_G = 17.82 * 10^-2 * kg; %(L/kg) Glucose distribution volume [Havorka et al 2002]
R_thr = 10.07; %(mmol/L)
R_cl = 1.05 * 10^-2; %(1/min)

S_IT = 11.9 * 10^-4; %(L/min/mU) Insulin sensitivity of distribution/transport [Havorka et al 2002]
S_ID = 6.64 * 10^-4; %(L/min/mU) Insulin sensitivity of disposal [Havorka et al 2002]
S_IE = 116 * 10^-4; %(L/min/mU) Insulin sensitivity of EGP [Havorka et al 2002]
k_b1 = .14 * 10^-2;
k_b2 = 13.77 * 10^-2;
k_b3 = 2.10 * 10^-2;
k_a1 = S_IT * k_b1; %(1/min) Deactivation rate [Havorka et al 2002]
k_a2 = S_ID * k_b2; %(1/min) Deactivation rate [Havorka et al 2002]
k_a3 = S_IE * k_b3; %(1/min) Deactivation rate [Havorka et al 2002]

k_a = 2.53 * 10^-2; %(1/min)
k_e = 14.0 * 10^-2; %(1/min) Insulin elimination from plasma [Havorka et al 1993]
V_I = 13.0 * 10^-2 * kg; %(L/kg) Insulin distribution volume [Havorka et al 1993]

Bio = .77; %(unitless) Carbohydrate Bioavailability [Livesey et al 1998]
t_max = 29; %(min) Time-to-maximum of CHO absorption [Livesey et al 1998]
U_G_ceil = .0275 * kg; %(mmol/kg/min) Maximum gut glucose flux drawn from U [.02, .035]

k_a_int = 11.03 * 10^-2;

t_maxI = 55; %(min) Time-to-maximum of absorption of subcutaneously injected short-acting insulin [Howey et al 1994], [Rave et al 1999]
%}

%% Havorka Patient 6 Specific Parameters
%{
kg = 107; %(kg) Patient Weight
needs = 41; %(U/day) patient insulin needs

Fs_01 = 10.28 * 10^-3 * kg; %(mmol/kg/min) Non-insulin-dependent glucose flux [Havorka et al 2002]
EGP_0 = 10.49 * 10^-3 * kg; %(mmol/kg/min) EGP extrapolated to zero insulin concentration [Havorka et al 2002]
k_12 = 3.07 * 10^-2; %(1/min)  Transfer rate [Havorka et al 2002]
V_G = 14.57 * 10^-2 * kg; %(L/kg) Glucose distribution volume [Havorka et al 2002]
R_thr = 7.75; %(mmol/L)
R_cl = 9.65 * 10^-2; %(1/min)

S_IT = 12.4 * 10^-4; %(L/min/mU) Insulin sensitivity of distribution/transport [Havorka et al 2002]
S_ID = 1.53 * 10^-4; %(L/min/mU) Insulin sensitivity of disposal [Havorka et al 2002]
S_IE = 114 * 10^-4; %(L/min/mU) Insulin sensitivity of EGP [Havorka et al 2002]
k_b1 = .07 * 10^-2;
k_b2 = 3.69 * 10^-2;
k_b3 = 3.39 * 10^-2;
k_a1 = S_IT * k_b1; %(1/min) Deactivation rate [Havorka et al 2002]
k_a2 = S_ID * k_b2; %(1/min) Deactivation rate [Havorka et al 2002]
k_a3 = S_IE * k_b3; %(1/min) Deactivation rate [Havorka et al 2002]

k_a = 2.57 * 10^-2; %(1/min)
k_e = 17.7 * 10^-2; %(1/min) Insulin elimination from plasma [Havorka et al 1993]
V_I = 11.5 * 10^-2 * kg; %(L/kg) Insulin distribution volume [Havorka et al 1993]

Bio = .72; %(unitless) Carbohydrate Bioavailability [Livesey et al 1998]
t_max = 26; %(min) Time-to-maximum of CHO absorption [Livesey et al 1998]
U_G_ceil = .0275 * kg; %(mmol/kg/min) Maximum gut glucose flux drawn from U [.02, .035]

k_a_int = 9.74 * 10^-2;

t_maxI = 55; %(min) Time-to-maximum of absorption of subcutaneously injected short-acting insulin [Howey et al 1994], [Rave et al 1999]
%}

%% Havorka Patient 8 Specific Parameters
%{
kg = 76; %(kg) Patient Weight
needs = 16; %(U/day) patient insulin needs

Fs_01 = 12.98 * 10^-3 * kg; %(mmol/kg/min) Non-insulin-dependent glucose flux [Havorka et al 2002]
EGP_0 = 14.07 * 10^-3 * kg; %(mmol/kg/min) EGP extrapolated to zero insulin concentration [Havorka et al 2002]
k_12 = 2.93 * 10^-2; %(1/min)  Transfer rate [Havorka et al 2002]
V_G = 18.03 * 10^-2 * kg; %(L/kg) Glucose distribution volume [Havorka et al 2002]
R_thr = 7.75; %(mmol/L)
R_cl = 1.15 * 10^-2; %(1/min)

S_IT = 30.0 * 10^-4; %(L/min/mU) Insulin sensitivity of distribution/transport [Havorka et al 2002]
S_ID = 1.70 * 10^-4; %(L/min/mU) Insulin sensitivity of disposal [Havorka et al 2002]
S_IE = 219 * 10^-4; %(L/min/mU) Insulin sensitivity of EGP [Havorka et al 2002]
k_b1 = .32 * 10^-2;
k_b2 = 21.95 * 10^-2;
k_b3 = 3.23 * 10^-2;
k_a1 = S_IT * k_b1; %(1/min) Deactivation rate [Havorka et al 2002]
k_a2 = S_ID * k_b2; %(1/min) Deactivation rate [Havorka et al 2002]
k_a3 = S_IE * k_b3; %(1/min) Deactivation rate [Havorka et al 2002]

k_a = 2.89 * 10^-2; %(1/min)
k_e = 12.0 * 10^-2; %(1/min) Insulin elimination from plasma [Havorka et al 1993]
V_I = 12.8 * 10^-2 * kg; %(L/kg) Insulin distribution volume [Havorka et al 1993]

Bio = .78; %(unitless) Carbohydrate Bioavailability [Livesey et al 1998]
t_max = 42; %(min) Time-to-maximum of CHO absorption [Livesey et al 1998]
U_G_ceil = .0275 * kg; %(mmol/kg/min) Maximum gut glucose flux drawn from U [.02, .035]

k_a_int = 6.89 * 10^-2;

t_maxI = 55; %(min) Time-to-maximum of absorption of subcutaneously injected short-acting insulin [Howey et al 1994], [Rave et al 1999]
%}

%% Havorka Patient 10 Specific Parameters
%{
kg = 101; %(kg) Patient Weight
needs = 32; %(U/day) patient insulin needs

Fs_01 = 3.96 * 10^-3 * kg; %(mmol/kg/min) Non-insulin-dependent glucose flux [Havorka et al 2002]
EGP_0 = 4.35 * 10^-3 * kg; %(mmol/kg/min) EGP extrapolated to zero insulin concentration [Havorka et al 2002]
k_12 = 5.37 * 10^-2; %(1/min)  Transfer rate [Havorka et al 2002]
V_G = 13.17 * 10^-2 * kg; %(L/kg) Glucose distribution volume [Havorka et al 2002]
R_thr = 7.55; %(mmol/L)
R_cl = 1.11 * 10^-2; %(1/min)

S_IT = 5.4 * 10^-4; %(L/min/mU) Insulin sensitivity of distribution/transport [Havorka et al 2002]
S_ID = 12.00 * 10^-4; %(L/min/mU) Insulin sensitivity of disposal [Havorka et al 2002]
S_IE = 53 * 10^-4; %(L/min/mU) Insulin sensitivity of EGP [Havorka et al 2002]
k_b1 = .48 * 10^-2;
k_b2 = 4.42 * 10^-2;
k_b3 = 1.66 * 10^-2;
k_a1 = S_IT * k_b1; %(1/min) Deactivation rate [Havorka et al 2002]
k_a2 = S_ID * k_b2; %(1/min) Deactivation rate [Havorka et al 2002]
k_a3 = S_IE * k_b3; %(1/min) Deactivation rate [Havorka et al 2002]

k_a = 2.44 * 10^-2; %(1/min)
k_e = 16.7 * 10^-2; %(1/min) Insulin elimination from plasma [Havorka et al 1993]
V_I = 12.9 * 10^-2 * kg; %(L/kg) Insulin distribution volume [Havorka et al 1993]

Bio = .71; %(unitless) Carbohydrate Bioavailability [Livesey et al 1998]
t_max = 52; %(min) Time-to-maximum of CHO absorption [Livesey et al 1998]
U_G_ceil = .0275 * kg; %(mmol/kg/min) Maximum gut glucose flux drawn from U [.02, .035]

k_a_int = 8.98 * 10^-2;

t_maxI = 55; %(min) Time-to-maximum of absorption of subcutaneously injected short-acting insulin [Howey et al 1994], [Rave et al 1999]
%}

%% Havorka Parameters from Normal Distribution
%{
kg = 90; %(kg) Patient Weight
needs = 39; %(U/day) patient insulin needs

% Generate Random Distribution
Fs_01 = 11.1 * 10^-3 * kg; %(mmol/kg/min) Non-insulin-dependent glucose flux [Havorka et al 2002]
EGP_0 = lognrnd(log(.0169), .86) * kg; %(mmol/kg/min) EGP extrapolated to zero insulin concentration [Havorka et al 2002]
k_12 = max(min(lognrnd(-2.813, .43), .2), .01); %(1/min)  Transfer rate [Havorka et al 2002]
V_G = max(min(lognrnd(-1.897, .23), .25), .09) * kg; %(L/kg) Glucose distribution volume [Havorka et al 2002]
R_thr = max(min(normrnd(9, 1.5), 15.0), 7.5); %(mmol/L)
R_cl = max(min(normrnd(.01, .025), .03), .003); %(1/min)

S_IT = 11.0 * 10^-4; %(L/min/mU) Insulin sensitivity of distribution/transport [Havorka et al 2002]
S_ID = 1.58 * 10^-4; %(L/min/mU) Insulin sensitivity of disposal [Havorka et al 2002]
S_IE = 73 * 10^-4; %(L/min/mU) Insulin sensitivity of EGP [Havorka et al 2002]
k_b1 = max(min(lognrnd(-5.684, 1.00), .05), .0002);
k_b2 = max(min(lognrnd(-2.882, .75), .4), .005);
k_b3 = max(min(lognrnd(-3.730, .75), .1), .003);
k_a1 = S_IT * k_b1; %(1/min) Deactivation rate [Havorka et al 2002]
k_a2 = S_ID * k_b2; %(1/min) Deactivation rate [Havorka et al 2002]
k_a3 = S_IE * k_b3; %(1/min) Deactivation rate [Havorka et al 2002]

out = [k_b1 k_b2 k_b3 k_a1 k_a2 k_a3]

k_a = max(min(normrnd(.018, .0045), .06), .005); %(1/min)
k_e = max(min(normrnd(.14, .0345), .3), .05); %(1/min) Insulin elimination from plasma [Havorka et al 1993]
V_I = max(min(normrnd(.12, .012), .18), .08) * kg; %(L/kg) Insulin distribution volume [Havorka et al 1993]

Bio = min(sqrt(120*12)*rand() + 50, 87) / 100; %(unitless) Carbohydrate Bioavailability [Livesey et al 1998] drawn from U(70,120)
t_max = 1 / max(min(lognrnd(-3.689, 0.252), .035), .02); %(min) Time-to-maximum of CHO absorption [Livesey et al 1998]
U_G_ceil = ((.035-.02)*rand() + .02) * kg; %(mmol/kg/min) Maximum gut glucose flux drawn from U [.02, .035]

k_a_int = lognrnd(-2.372, 1.092);

t_maxI = 55; %(min) Time-to-maximum of absorption of subcutaneously injected short-acting insulin [Howey et al 1994], [Rave et al 1999]
save('random_vars.mat')
%}

%% Initial Conditions
% G_0 = 25;%((30-2)*rand() + 2); %(mmol/l) Uniform distribution for glucose concentration
Q_1_0 = G_0 * V_G;
%I_0_p = Fs_01 / (2*S_IT*(Q_1_0 + V_G)) * (-1 + sqrt(1 - (4*S_IT*k_12*(Q_1_0 * V_G)) / (Fs_01 * S_ID)));
%I_0_m = Fs_01 / (2*S_IT*(Q_1_0 + V_G)) * (-1 - sqrt(1 - (4*S_IT*k_12*(Q_1_0 * V_G)) / (Fs_01 * S_ID)));
I_0 = 1.15;

x_1_0 = S_IT * I_0;
x_2_0 = S_ID * I_0;
x_3_0 = S_IE * I_0;

Q_2_0 = x_1_0 / (k_12 + x_2_0) * Q_1_0;

G_1_0 = 0;
G_2_0 = 0;

S_1_0 = 0;
S_2_0 = 0;

% %% Brute Force MPC Optimization
% 
% dose_start = 500;
% control_steps = 100;
% step_length = 1; %(mins)
% horiz_num = 1;
% pred_horiz = 400 * horiz_num; %(mins)
% 
% % Model Information
% paramNameValStruct.SimulationMode = 'normal';
% paramNameValStruct.AbsTol         = '1e-5';
% paramNameValStruct.SaveState      = 'on';
% paramNameValStruct.StateSaveName  = 'xout';
% paramNameValStruct.SaveOutput     = 'on';
% paramNameValStruct.OutputSaveName = 'yout';
% paramNameValStruct.SaveFormat = 'Dataset';
% paramNameValStruct.StartTime = '0';
% 
% open_system('Havorka_2004_distribution');
% open_system('Havorka_2004_distribution_final');
% 
% insulin_sweep = 0:.5:10; %(mU) Insulin dose sweeping range
% %get_param('Havorka_2004_distribution/dose', 'ObjectParameters')
% %get_param('Havorka_2004_distribution/dose', 'DialogParameters')
% insulin_dose = zeros(control_steps);
% 
% glucose_final = zeros(length(insulin_sweep));
% glucose_target = 0;
% 
% %figure;
% 
% for n = 0 : (control_steps - 1)
%     
%     if n == 0
%         set_param('Havorka_2004_distribution_final/dose', 'Value', 'insulin_dose(n+1)');
%         simFinal = sim('Havorka_2004_distribution_final', 'AbsTol', '1e-5', ...
%             'StartTime', '0', 'StopTime', 'dose_start',  ...
%             'SaveOutput', 'on', 'OutputSaveName', 'yFinal', 'SaveState', 'on', 'StateSaveName', 'xFinal', ...
%             'SaveFinalState', 'on', 'FinalStateName', 'stateFinal', 'SaveFormat', 'Structure');
%         
%         stateFinal = simFinal.get('stateFinal');
%         yFinal = simFinal.get('yFinal');
% %         glucose = (yFinal.get('glucose').Values);
%         glucose = yFinal.signals(1).values;
%         glucose_target = glucose(end);
%         
%     else
%         set_param('Havorka_2004_distribution_final/dose', 'Value', 'insulin_dose(n-1)');
%         simFinal = sim('Havorka_2004_distribution_final', 'AbsTol', '1e-5', ...
%             'StartTime', 'dose_start + (n-1)*step_length', 'StopTime', 'dose_start + (n)*step_length', ...
%             'SaveOutput', 'on', 'OutputSaveName', 'yfinal', 'SaveState', 'on', 'StateSaveName', 'xfinal', ...
%             'SaveFinalState', 'on', 'FinalStateName', 'stateFinal', 'SaveFormat', 'Structure');
%         stateFinal = simFinal.get('stateFinal');
%         yFinal = simFinal.get('yFinal');
% %         glucose = (yFinal.get('glucose').Values);
%         glucose = yFinal.signals(1).values;
%         glucose_target = glucose(end);
%         
%      end
%     
%     for i = 1:length(insulin_sweep)
%         % Run the Tested Insulin Dose
%         set_param('Havorka_2004_distribution/dose', 'Value', 'insulin_sweep(i)');
%         simOut = sim('Havorka_2004_distribution', 'AbsTol', '1e-5', ...
%             'StartTime', 'dose_start + n*step_length', 'StopTime', 'dose_start + (n+1)*step_length', ...
%             'SaveOutput', 'on', 'OutputSaveName', 'yout', 'SaveState', 'on', 'StateSaveName', 'xout', ...
%             'SaveFinalState', 'on', 'FinalStateName', 'stateOut', ...
%             'LoadInitialState', 'on', 'InitialState', 'stateFinal', 'SaveFormat', 'Structure');
%         
%         % Obtain Model Outputs
%         outputs = simOut.get('yout');
% %         glucose =(outputs.get('glucose').Values);
% %         insulin =(outputs.get('insulin').Values);
% %         controller = (outputs.get('controller').Values);
%         glucose = outputs.signals(1).values;
%         insulin = outputs.signals(2).values;
%         controller = outputs.signals(3).values;
%         
%         glucose_final(i) = glucose(end);
%         
%         % Decide on Optimal Insulin Dose
%         if (i > 1) & (glucose_final(i) >= glucose_target) & (glucose_final(i) <= glucose_final(i-1))
%             insulin_dose = insulin_sweep(i);
%         
%         else
%             insulin_dose = insulin_sweep(i);
%         end
%         
%         % Plot Model Outputs
%         %{
%         subplot(3,1,1)
%         plot1 = plot(glucose);
%         set(plot1, 'color', [0 1 0]);
%         axis([dose_start dose_start+pred_horiz 10 30]);
%         title('Glucose');
%         hold on;
% 
%         subplot(3,1,2)
%         plot2 = plot(insulin);
%         set(plot2, 'color', [0 0 1]);
%         axis([dose_start dose_start+pred_horiz 0 10]);
%         title('Insulin');
%         hold on;
% 
%         subplot(3,1,3)
%         plot3 = plot(controller);
%         set(plot3, 'color', [1 0 0]);
%         axis([dose_start dose_start+pred_horiz 0 6]);
%         title('Controller Insulin Dose');
%         hold on;
%         %}
%         
%         disp(i)
%     end
% end