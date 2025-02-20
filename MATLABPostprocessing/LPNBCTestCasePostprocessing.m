%% LPN Boundary Condition Test Case Postprocessing Script

%{
Author:
    Muhammad Ahmad Raza, University College Dublin, Ireland.

Cite as:
    Raza, M. A.: Implementation of Lumped Parameter Network Boundary Conditions for the
    Patient-Specific CFD Simulations of Coronary Arteries in OpenFOAM. In Proceedings of CFD with
    OpenSource Software, 2024, Edited by Nilsson. H., http://dx.doi.org/10.17196/OS_CFD#YEAR_2024
%}

%% Clear All

clc;
clear all;
close all;
%% Parameters Values

N = 8; %Number of cycles

%Resistive model parameters
R_Rd = 17690000;
R_Pd = 9.7363e+03;

%2-Element Windkessel model parameters
WK2_Rd = 141520000;
WK2_C = 8.3333e-09;
WK2_Pd = 0;

%3-Element Windkessel model parameters
WK3_Rp = 13997000;
WK3_Rd = 141520000;
WK3_C = 1.0000e-08;
WK3_Pd = 0;

%4-Element (Series) Windkessel model parameters
WK4s_Rp = 13997000;
WK4s_Rd = 141520000;
WK4s_C = 2.5000e-08;
WK4s_L = 1.7995e+04;
WK4s_Pd = 0;

%4-Element (Parallel) Windkessel model parameters
WK4p_Rp = 1.3997e+10;
WK4p_Rd = 141520000;
WK4p_C = 1.0000e-08;
WK4p_L = 1.7995e+03;
WK4p_Pd = 0;

%Coronary LPN Model parameters
Ra = 1.0179e+10; %Arterial resistance [kg m^-4 s^-1]
Ram = 1.6541e+10; %Micro-arterial resistance [kg m^-4 s^-1]
Rv = 5.0896e+09; %Veinous resistance [kg m^-4 s^-1]
Rvm = 0; %Micro-veinous resistance [kg m^-4 s^-1]
Ca = 8.6482e-12; %Arterial compliance [m^4 s^2 kg^-1]
Cim = 6.9972e-11; %Intramyocardial compliance [m^4 s^2 kg^-1]
PimScaling = 1.5; %Intramyocardial pressure scaling (1.5 for LCA, 0.5 for RCA)
Pv = 0; %Distal pressure [Pa]

%% Read and prepare Aortic Flow Rate Data

% Open the flow rate data file
AorticFlowfilename = '../LPNBCTestCase/DataFiles/AorticInletFlowRate'; % Replace with your actual file name
AorticFlowfileID = fopen(AorticFlowfilename, 'r');

% Read the entire content as a string
AorticFlowrawData = fscanf(AorticFlowfileID, '%c');

% Close the file
fclose(AorticFlowfileID);

% Remove outer parentheses
AorticFlowcleanedData = strrep(AorticFlowrawData, '(', '');
AorticFlowcleanedData = strrep(AorticFlowcleanedData, ')', '');

% Convert to numeric array
AorticFlowData = sscanf(AorticFlowcleanedData, '%f %f', [2, Inf])';

%%
% Repeat the Aortic Flow Rate data for multiple cycles
cycleTime = max(AorticFlowData(:,1)); % Duration of one cycle
AorticFlowTime = []; % Initialize extended time array
AorticFlowRate = []; % Initialize extended flow rate array

for i = 0:(N-1)
    % Offset time for the i-th cycle and avoid duplication at overlap
    newAorticFlowTime = AorticFlowData(:,1) + i * cycleTime;
    newAorticFlowRate = AorticFlowData(:,2);
    if i > 0
        newAorticFlowTime(1) = []; % Remove the first point to avoid duplicate
        newAorticFlowRate(1) = []; % Remove the corresponding flow rate value
    end
    AorticFlowTime = [AorticFlowTime;  newAorticFlowTime];
    AorticFlowRate = [AorticFlowRate;  newAorticFlowRate];
end

%% Read and prepare Coronary Flow Rate Data

% Open the flow rate data file
CoronaryFlowfilename = '../LPNBCTestCase/DataFiles/CoronaryInletFlowRate'; % Replace with your actual file name
CoronaryFlowfileID = fopen(CoronaryFlowfilename, 'r');

% Read the entire content as a string
CoronaryFlowrawData = fscanf(CoronaryFlowfileID, '%c');

% Close the file
fclose(CoronaryFlowfileID);

% Remove outer parentheses
CoronaryFlowcleanedData = strrep(CoronaryFlowrawData, '(', '');
CoronaryFlowcleanedData = strrep(CoronaryFlowcleanedData, ')', '');

% Convert to numeric array
CoronaryFlowData = sscanf(CoronaryFlowcleanedData, '%f %f', [2, Inf])';

%%
% Repeat the Coronary Flow Rate data for multiple cycles
cycleTime = max(CoronaryFlowData(:,1)); % Duration of one cycle
CoronaryFlowTime = []; % Initialize extended time array
CoronaryFlowRate = []; % Initialize extended flow rate array

for i = 0:(N-1)
    % Offset time for the i-th cycle and avoid duplication at overlap
    newCoronaryFlowTime = CoronaryFlowData(:,1) + i * cycleTime;
    newCoronaryFlowRate = CoronaryFlowData(:,2);
    if i > 0
        newCoronaryFlowTime(1) = []; % Remove the first point to avoid duplicate
        newCoronaryFlowRate(1) = []; % Remove the corresponding flow rate value
    end
    CoronaryFlowTime = [CoronaryFlowTime;  newCoronaryFlowTime];
    CoronaryFlowRate = [CoronaryFlowRate;  newCoronaryFlowRate];
end

%% Read and prepare Pim Data File

% Open the Pim data file
Pimfilename = '../LPNBCTestCase/DataFiles/PimData'; % Replace with your actual file name
PimfileID = fopen(Pimfilename, 'r');

% Read the entire content as a string
PimrawData = fscanf(PimfileID, '%c');

% Close the file
fclose(PimfileID);

% Remove outer parentheses
PimcleanedData = strrep(PimrawData, '(', '');
PimcleanedData = strrep(PimcleanedData, ')', '');

% Convert to numeric array
PimData = sscanf(PimcleanedData, '%f %f', [2, Inf])';

%%
% Repeat the Coronary Flow Rate data for multiple cycles
PimcycleTime = max(PimData(:,1)); % Duration of one cycle
PimTime = []; % Initialize extended time array
Pim = []; % Initialize extended flow rate array

for i = 0:(N-1)
    % Offset time for the i-th cycle and avoid duplication at overlap
    newPimTime = PimData(:,1) + i * PimcycleTime;
    newPim = PimData(:,2);
    if i > 0
        newPimTime(1) = []; % Remove the first point to avoid duplicate
        newPim(1) = []; % Remove the corresponding flow rate value
    end
    PimTime = [PimTime;  newPimTime];
    Pim = [Pim;  newPim];
end

%% Solve Windkessel Models

[R_tSol, R_PSol] = Resistive(R_Rd, R_Pd, AorticFlowRate, AorticFlowTime);

[WK2_tSol, WK2_PSol]= WK2(WK2_Rd, WK2_C, WK2_Pd, AorticFlowRate, AorticFlowTime);

[WK3_tSol, WK3_PSol] = WK3(WK3_Rp, WK3_Rd, WK3_C, WK3_Pd, AorticFlowRate, AorticFlowTime);

[WK4s_tSol, WK4s_PSol] = WK4Series(WK4s_Rp, WK4s_Rd, WK4s_C, WK4s_L, WK4s_Pd, AorticFlowRate, AorticFlowTime);

[WK4p_tSol, WK4p_PSol] = WK4Parallel(WK4p_Rp, WK4p_Rd, WK4p_C, WK4p_L, WK4p_Pd, AorticFlowRate, AorticFlowTime);

%% Solve Coronary LPN Model

[Coronary_tSol, Coronary_PSol] = CoronaryLPN(Ra, Ram, Rv, Rvm, Ca, Cim, Pv, PimScaling, Pim, PimTime, CoronaryFlowRate, CoronaryFlowTime);

%% Read simulation data

R_Qout = table2array(readtable('../LPNBCTestCase_Resistive_/postProcessing/flowRateOutlet/0/faceSource.dat'));
R_Pout = table2array(readtable('../LPNBCTestCase_Resistive_/postProcessing/pAverageOutlet/0/faceSource.dat'));

WK2_Qout1 = table2array(readtable('../LPNBCTestCase_WK2_firstOrder/postProcessing/flowRateOutlet/0/faceSource.dat'));
WK2_Qout2 = table2array(readtable('../LPNBCTestCase_WK2_secondOrder/postProcessing/flowRateOutlet/0/faceSource.dat'));
WK2_Pout1 = table2array(readtable('../LPNBCTestCase_WK2_firstOrder/postProcessing/pAverageOutlet/0/faceSource.dat'));
WK2_Pout2 = table2array(readtable('../LPNBCTestCase_WK2_secondOrder/postProcessing/pAverageOutlet/0/faceSource.dat'));

WK3_Qout1 = table2array(readtable('../LPNBCTestCase_WK3_firstOrder/postProcessing/flowRateOutlet/0/faceSource.dat'));
WK3_Qout2 = table2array(readtable('../LPNBCTestCase_WK3_secondOrder/postProcessing/flowRateOutlet/0/faceSource.dat'));
WK3_Pout1 = table2array(readtable('../LPNBCTestCase_WK3_firstOrder/postProcessing/pAverageOutlet/0/faceSource.dat'));
WK3_Pout2 = table2array(readtable('../LPNBCTestCase_WK3_secondOrder/postProcessing/pAverageOutlet/0/faceSource.dat'));

WK4s_Qout1 = table2array(readtable('../LPNBCTestCase_WK4Series_firstOrder/postProcessing/flowRateOutlet/0/faceSource.dat'));
WK4s_Qout2 = table2array(readtable('../LPNBCTestCase_WK4Series_secondOrder/postProcessing/flowRateOutlet/0/faceSource.dat'));
WK4s_Pout1 = table2array(readtable('../LPNBCTestCase_WK4Series_firstOrder/postProcessing/pAverageOutlet/0/faceSource.dat'));
WK4s_Pout2 = table2array(readtable('../LPNBCTestCase_WK4Series_secondOrder/postProcessing/pAverageOutlet/0/faceSource.dat'));

WK4p_Qout1 = table2array(readtable('../LPNBCTestCase_WK4Parallel_firstOrder/postProcessing/flowRateOutlet/0/faceSource.dat'));
WK4p_Qout2 = table2array(readtable('../LPNBCTestCase_WK4Parallel_secondOrder/postProcessing/flowRateOutlet/0/faceSource.dat'));
WK4p_Pout1 = table2array(readtable('../LPNBCTestCase_WK4Parallel_firstOrder/postProcessing/pAverageOutlet/0/faceSource.dat'));
WK4p_Pout2 = table2array(readtable('../LPNBCTestCase_WK4Parallel_secondOrder/postProcessing/pAverageOutlet/0/faceSource.dat'));

Coronary_Qout1 = table2array(readtable('../LPNBCTestCase_Coronary_firstOrder/postProcessing/flowRateOutlet/0/faceSource.dat'));
Coronary_Qout2 = table2array(readtable('../LPNBCTestCase_Coronary_secondOrder/postProcessing/flowRateOutlet/0/faceSource.dat'));
Coronary_Pout1 = table2array(readtable('../LPNBCTestCase_Coronary_firstOrder/postProcessing/pAverageOutlet/0/faceSource.dat'));
Coronary_Pout2 = table2array(readtable('../LPNBCTestCase_Coronary_secondOrder/postProcessing/pAverageOutlet/0/faceSource.dat'));

%% Plot results

% Plot Resistive Model results
figure('defaultAxesFontSize', 12, 'defaultLineLineWidth', 1, ...
       'Units', 'inches', 'Position', [1, 1, 9, 6]);

% Define layout for the subplots
tiledlayout(3, 1); % 3 rows, 1 column (upper subplot takes 2 rows, bottom takes 1 row)

% Upper subplot for Pressure
nexttile([2 1]); % Occupy 2 rows
hold on;

plot(R_tSol, R_PSol / 133.33, '-k', 'LineWidth', 1.5); % Pressure plot
plot(R_Pout(:,1), R_Pout(:,2) /133.33, '--r', 'LineWidth', 1.5);
xlabel('$t \ (s)$', 'Interpreter', 'latex');
ylabel('$P(t) \ (mmHg)$', 'Interpreter', 'latex');
%title('Pressure Solution', 'Interpreter', 'latex');
legend('MATLAB', 'Resistive BC', 'Interpreter', 'latex', 'location', 'southeast', 'Box', 'off');
ax = gca;
ax.YColor = 'k'; % Set left axis color to blue
ax.TickLabelInterpreter = 'latex'; % Use LaTeX for tick labels
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
grid on; grid minor; box on;
% Adjust figure for better visibility
set(gca, 'FontSize', 12);
xlim([0 8]); ylim([70 120]);
hold off;

% Lower subplot for Flow Rate
nexttile; % Occupy 1 row
hold on;

plot(AorticFlowTime, AorticFlowRate * 1e6, '-k', 'LineWidth', 1.5); % Flow rate plot
plot(R_Qout(:,1), R_Qout(:,2) * 1e6, '--r', 'LineWidth', 1.5);
xlabel('$t \ (s)$', 'Interpreter', 'latex');
ylabel('$Q(t) \ (mL/s)$', 'Interpreter', 'latex');
%title('Flow Rate Solution', 'Interpreter', 'latex');
legend('Inlet','Resistive BC', 'Interpreter','latex', 'location', 'southeast', 'Box', 'off')
ax = gca;
ax.YColor = 'k'; % Set left axis color to red
ax.TickLabelInterpreter = 'latex'; % Use LaTeX for tick labels
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
grid on; grid minor; box on;
% Adjust figure for better visibility
set(gca, 'FontSize', 12);
xlim([0 8]); ylim([-50 350]);
% Save the figure as an EPS file
print(gcf, 'R_Pressure', '-depsc', '-r300'); % '-depsc' for color EPS, '-r300' for resolution
hold off; 


% Plot 2-Element Windkessel Model results
figure('defaultAxesFontSize', 12, 'defaultLineLineWidth', 1, ...
       'Units', 'inches', 'Position', [1, 1, 9, 6]);

% Define layout for the subplots
tiledlayout(3, 1); % 3 rows, 1 column (upper subplot takes 2 rows, bottom takes 1 row)

% Upper subplot for Pressure
nexttile([2 1]); % Occupy 2 rows
hold on;

plot(WK2_tSol, WK2_PSol / 133.33, '-k', 'LineWidth', 1.5); % Pressure plot
plot(WK2_Pout1(:,1), WK2_Pout1(:,2) /133.33, '--r', 'LineWidth', 1.5);
plot(WK2_Pout2(:,1), WK2_Pout2(:,2) /133.33, '--b', 'LineWidth', 1.5);
xlabel('$t \ (s)$', 'Interpreter', 'latex');
ylabel('$P(t) \ (mmHg)$', 'Interpreter', 'latex');
%title('Pressure Solution', 'Interpreter', 'latex');
legend('MATLAB', '$O(\Delta t)$ WK2 BC', '$O(\Delta t^2)$ WK2 BC', 'Interpreter', 'latex', 'location', 'southeast', 'Box', 'off');
ax = gca;
ax.YColor = 'k'; % Set left axis color to blue
ax.TickLabelInterpreter = 'latex'; % Use LaTeX for tick labels
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
grid on; grid minor; box on;
% Adjust figure for better visibility
set(gca, 'FontSize', 12);
xlim([0 8]); ylim([0 130]);
hold off;

% Lower subplot for Flow Rate
nexttile; % Occupy 1 row
hold on;

plot(AorticFlowTime, AorticFlowRate * 1e6, '-k', 'LineWidth', 1.5); % Flow rate plot
plot(WK2_Qout1(:,1), WK2_Qout1(:,2) * 1e6, '--r', 'LineWidth', 1.5);
plot(WK2_Qout2(:,1), WK2_Qout2(:,2) * 1e6, '--b', 'LineWidth', 1.5);
xlabel('$t \ (s)$', 'Interpreter', 'latex');
ylabel('$Q(t) \ (mL/s)$', 'Interpreter', 'latex');
%title('Flow Rate Solution', 'Interpreter', 'latex');
legend('Inlet','$O(\Delta t)$ WK2 BC', '$O(\Delta t^2)$ WK2 BC', 'Interpreter','latex', 'location', 'southeast', 'Box', 'off')
ax = gca;
ax.YColor = 'k'; % Set left axis color to red
ax.TickLabelInterpreter = 'latex'; % Use LaTeX for tick labels
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
grid on; grid minor; box on;
% Adjust figure for better visibility
set(gca, 'FontSize', 12);
xlim([0 8]); ylim([-50 350]);
% Save the figure as an EPS file
print(gcf, 'WK2_Pressure', '-depsc', '-r300'); % '-depsc' for color EPS, '-r300' for resolution
hold off;

% Plot 3-Element Windkessel Model results
figure('defaultAxesFontSize', 12, 'defaultLineLineWidth', 1, ...
       'Units', 'inches', 'Position', [1, 1, 9, 6]);

% Define layout for the subplots
tiledlayout(3, 1); % 3 rows, 1 column (upper subplot takes 2 rows, bottom takes 1 row)

% Upper subplot for Pressure
nexttile([2 1]); % Occupy 2 rows
hold on;

plot(WK3_tSol, WK3_PSol / 133.33, '-k', 'LineWidth', 1.5); % Pressure plot
plot(WK3_Pout1(:,1), WK3_Pout1(:,2) /133.33, '--r', 'LineWidth', 1.5);
plot(WK3_Pout2(:,1), WK3_Pout2(:,2) /133.33, '--b', 'LineWidth', 1.5);
xlabel('$t \ (s)$', 'Interpreter', 'latex');
ylabel('$P(t) \ (mmHg)$', 'Interpreter', 'latex');
%title('Pressure Solution', 'Interpreter', 'latex');
legend('MATLAB', '$O(\Delta t)$ WK3 BC', '$O(\Delta t^2)$ WK3 BC', 'Interpreter', 'latex', 'location', 'southeast', 'Box', 'off');
ax = gca;
ax.YColor = 'k'; % Set left axis color to blue
ax.TickLabelInterpreter = 'latex'; % Use LaTeX for tick labels
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
grid on; grid minor; box on;
% Adjust figure for better visibility
set(gca, 'FontSize', 12);
xlim([0 8]); ylim([0 140]);
hold off;

% Lower subplot for Flow Rate
nexttile; % Occupy 1 row
hold on;

plot(AorticFlowTime, AorticFlowRate * 1e6, '-k', 'LineWidth', 1.5); % Flow rate plot
plot(WK3_Qout1(:,1), WK3_Qout1(:,2) * 1e6, '--r', 'LineWidth', 1.5);
plot(WK3_Qout2(:,1), WK3_Qout2(:,2) * 1e6, '--b', 'LineWidth', 1.5);
xlabel('$t \ (s)$', 'Interpreter', 'latex');
ylabel('$Q(t) \ (mL/s)$', 'Interpreter', 'latex');
%title('Flow Rate Solution', 'Interpreter', 'latex');
legend('Inlet','$O(\Delta t)$ WK3 BC', '$O(\Delta t^2)$ WK3 BC', 'Interpreter','latex', 'location', 'southeast', 'Box', 'off')
ax = gca;
ax.YColor = 'k'; % Set left axis color to red
ax.TickLabelInterpreter = 'latex'; % Use LaTeX for tick labels
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
grid on; grid minor; box on;
% Adjust figure for better visibility
set(gca, 'FontSize', 12);
xlim([0 8]); ylim([-50 350]);
% Save the figure as an EPS file
print(gcf, 'WK3_Pressure', '-depsc', '-r300'); % '-depsc' for color EPS, '-r300' for resolution
hold off;


% Plot 4-Element (Series) Windkessel Model results
figure('defaultAxesFontSize', 12, 'defaultLineLineWidth', 1, ...
       'Units', 'inches', 'Position', [1, 1, 9, 6]);

% Define layout for the subplots
tiledlayout(3, 1); % 3 rows, 1 column (upper subplot takes 2 rows, bottom takes 1 row)

% Upper subplot for Pressure
nexttile([2 1]); % Occupy 2 rows
hold on;

plot(WK4s_tSol, WK4s_PSol / 133.33, '-k', 'LineWidth', 1.5); % Pressure plot
plot(WK4s_Pout1(:,1), WK4s_Pout1(:,2) /133.33, '--r', 'LineWidth', 1.5);
plot(WK4s_Pout2(:,1), WK4s_Pout2(:,2) /133.33, '--b', 'LineWidth', 1.5);
xlabel('$t \ (s)$', 'Interpreter', 'latex');
ylabel('$P(t) \ (mmHg)$', 'Interpreter', 'latex');
%title('Pressure Solution', 'Interpreter', 'latex');
legend('MATLAB', '$O(\Delta t)$ WK4Series BC', '$O(\Delta t^2)$ WK4Series BC', 'Interpreter', 'latex', 'location', 'southeast', 'Box', 'off');
ax = gca;
ax.YColor = 'k'; % Set left axis color to blue
ax.TickLabelInterpreter = 'latex'; % Use LaTeX for tick labels
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
grid on; grid minor; box on;
% Adjust figure for better visibility
set(gca, 'FontSize', 12);
xlim([0 8]); ylim([0 130]);
hold off;

% Lower subplot for Flow Rate
nexttile; % Occupy 1 row
hold on;

plot(AorticFlowTime, AorticFlowRate * 1e6, '-k', 'LineWidth', 1.5); % Flow rate plot
plot(WK4s_Qout1(:,1), WK4s_Qout1(:,2) * 1e6, '--r', 'LineWidth', 1.5);
plot(WK4s_Qout2(:,1), WK4s_Qout2(:,2) * 1e6, '--b', 'LineWidth', 1.5);
xlabel('$t \ (s)$', 'Interpreter', 'latex');
ylabel('$Q(t) \ (mL/s)$', 'Interpreter', 'latex');
%title('Flow Rate Solution', 'Interpreter', 'latex');
legend('Inlet','$O(\Delta t)$ WK4Series BC', '$O(\Delta t^2)$ WK4Series BC', 'Interpreter','latex', 'location', 'southeast', 'Box', 'off')
ax = gca;
ax.YColor = 'k'; % Set left axis color to red
ax.TickLabelInterpreter = 'latex'; % Use LaTeX for tick labels
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
grid on; grid minor; box on;
% Adjust figure for better visibility
set(gca, 'FontSize', 12);
xlim([0 8]); ylim([-50 350]);
% Save the figure as an EPS file
print(gcf, 'WK4s_Pressure', '-depsc', '-r300'); % '-depsc' for color EPS, '-r300' for resolution
hold off;


% Plot 4-Element (Parallel) Windkessel Model results
figure('defaultAxesFontSize', 12, 'defaultLineLineWidth', 1, ...
       'Units', 'inches', 'Position', [1, 1, 9, 6]);

% Define layout for the subplots
tiledlayout(3, 1); % 3 rows, 1 column (upper subplot takes 2 rows, bottom takes 1 row)

% Upper subplot for Pressure
nexttile([2 1]); % Occupy 2 rows
hold on;

plot(WK4p_tSol, WK4p_PSol(:,1) / 133.33, '-k', 'LineWidth', 1.5); % Pressure plot
plot(WK4p_Pout1(:,1), WK4p_Pout1(:,2) /133.33, '--r', 'LineWidth', 1.5);
plot(WK4p_Pout2(:,1), WK4p_Pout2(:,2) /133.33, '--b', 'LineWidth', 1.5);
xlabel('$t \ (s)$', 'Interpreter', 'latex');
ylabel('$P(t) \ (mmHg)$', 'Interpreter', 'latex');
%title('Pressure Solution', 'Interpreter', 'latex');
legend('MATLAB', '$O(\Delta t)$ WK4Parallel BC', '$O(\Delta t^2)$ WK4Parallel BC', 'Interpreter', 'latex', 'location', 'southeast', 'Box', 'off');
ax = gca;
ax.YColor = 'k'; % Set left axis color to blue
ax.TickLabelInterpreter = 'latex'; % Use LaTeX for tick labels
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
grid on; grid minor; box on;
% Adjust figure for better visibility
set(gca, 'FontSize', 12);
xlim([0 8]); ylim([0 120]);
hold off;

% Lower subplot for Flow Rate
nexttile; % Occupy 1 row
hold on;

plot(AorticFlowTime, AorticFlowRate * 1e6, '-k', 'LineWidth', 1.5); % Flow rate plot
plot(WK4p_Qout1(:,1), WK4p_Qout1(:,2) * 1e6, '--r', 'LineWidth', 1.5);
plot(WK4p_Qout2(:,1), WK4p_Qout2(:,2) * 1e6, '--b', 'LineWidth', 1.5);
xlabel('$t \ (s)$', 'Interpreter', 'latex');
ylabel('$Q(t) \ (mL/s)$', 'Interpreter', 'latex');
%title('Flow Rate Solution', 'Interpreter', 'latex');
legend('Inlet','$O(\Delta t)$ WK4Parallel BC', '$O(\Delta t^2)$ WK4Parallel BC', 'Interpreter','latex', 'location', 'southeast', 'Box', 'off')
ax = gca;
ax.YColor = 'k'; % Set left axis color to red
ax.TickLabelInterpreter = 'latex'; % Use LaTeX for tick labels
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
grid on; grid minor; box on;
% Adjust figure for better visibility
set(gca, 'FontSize', 12);
xlim([0 8]); ylim([-50 350]);
% Save the figure as an EPS file
print(gcf, 'WK4p_Pressure', '-depsc', '-r300'); % '-depsc' for color EPS, '-r300' for resolution
hold off;

% Plot Coronary LPN Model results
figure('defaultAxesFontSize', 12, 'defaultLineLineWidth', 1, ...
       'Units', 'inches', 'Position', [1, 1, 9, 6]);

% Define layout for the subplots
tiledlayout(3, 1); % 3 rows, 1 column (upper subplot takes 2 rows, bottom takes 1 row)

% Upper subplot for Pressure
nexttile([2 1]); % Occupy 2 rows
hold on;

% Plotting on the left y-axis
yyaxis left;
plot(Coronary_tSol, Coronary_PSol(:,1) / 133.33, '-k', 'LineWidth', 1.5); % Pressure plot
plot(Coronary_Pout1(:,1), Coronary_Pout1(:,2) / 133.33, '--r', 'LineWidth', 1.5);
plot(Coronary_Pout2(:,1), Coronary_Pout2(:,2) / 133.33, '--b', 'LineWidth', 1.5);

ylabel('$P(t) \ (mmHg)$', 'Interpreter', 'latex'); % Left y-axis label
ax = gca;
ax.YColor = 'k'; % Set left axis color
ax.TickLabelInterpreter = 'latex'; % Use LaTeX for tick labels

% Plotting on the right y-axis
yyaxis right;
plot(PimTime, PimScaling * Pim / 133.33, ':', 'LineWidth', 1.5);
ylabel('$P_{im}(t) \ (mmHg)$', 'Interpreter', 'latex'); % Right y-axis label
ax = gca;
ax.YColor = [0.8500 0.3250 0.0980]; %'m'; % Set right axis color

xlabel('$t \ (s)$', 'Interpreter', 'latex'); % X-axis label
%title('Pressure Solution', 'Interpreter', 'latex');

legend('MATLAB', '$O(\Delta t)$ Coronary LPN BC', '$O(\Delta t^2)$ Coronary LPN BC', '$P_{im}(t)$', ...
       'Interpreter', 'latex', 'location', 'southeast', 'Box', 'off');

% Use LaTeX for tick labels and customize ticks
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
grid on; grid minor; box on;

% Adjust figure for better visibility
set(gca, 'FontSize', 12);
xlim([0 8]); 
yyaxis left; ylim([0 140]); % Adjust left y-axis limits
yyaxis right; ylim([0 200]); %ylim([min(PimScaling * Pim / 133.33), max(PimScaling * Pim / 133.33)]); % Adjust right y-axis limits

hold off;

% Lower subplot for Flow Rate
nexttile; % Occupy 1 row
hold on;

plot(CoronaryFlowTime, CoronaryFlowRate * 1e6, '-k', 'LineWidth', 1.5); % Flow rate plot
plot(Coronary_Qout1(:,1), Coronary_Qout1(:,2) * 1e6, '--r', 'LineWidth', 1.5);
plot(Coronary_Qout2(:,1), Coronary_Qout2(:,2) * 1e6, '--b', 'LineWidth', 1.5);
xlabel('$t \ (s)$', 'Interpreter', 'latex');
ylabel('$Q(t) \ (mL/s)$', 'Interpreter', 'latex');
%title('Flow Rate Solution', 'Interpreter', 'latex');
legend('Inlet','$O(\Delta t)$ Coronary LPN BC', '$O(\Delta t^2)$ Coronary LPN BC', 'Interpreter','latex', 'location', 'southeast', 'Box', 'off')
ax = gca;
ax.YColor = 'k'; % Set left axis color to red
ax.TickLabelInterpreter = 'latex'; % Use LaTeX for tick labels
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');
grid on; grid minor; box on;
% Adjust figure for better visibility
set(gca, 'FontSize', 12);
xlim([0 8]); ylim([0.15 0.7]);
% Save the figure as an EPS file
print(gcf, 'Coronary_Pressure', '-depsc', '-r300'); % '-depsc' for color EPS, '-r300' for resolution
hold off;

%% WK1 Function
function [tSol, PSol] = Resistive(Rd, Pd, FlowRate, FlowTime)

    % Q_fun represents the flow rate as a function of time
    Q_fun =@(t) interp1(FlowTime, FlowRate, t, 'linear', 'extrap');

    % Time span for the solution
    tSol = FlowTime;

    % Calculate P(t) using the WK1 equation
    PSol = Rd .* arrayfun(Q_fun, tSol) + Pd;

end


%% WK2 Function
function [tSol, PSol]= WK2(Rd, C, Pd, FlowRate, FlowTime)

    Q_fun =@(t) interp1(FlowTime, FlowRate, t, 'linear', 'extrap');

    % Define the differential equation
    % P'(t) + (1/Rd*C) * (P(t)-Pd) = Q(t)/C
    % Rearrange as: P'(t) = Q(t)/C - (1/Rd*C) * (P(t)-Pd)
    
    % Solve using MATLAB's ode45
    odefun = @(t, P) Q_fun(t) / C - (1 / (Rd * C)) * (P-Pd);
    
    % Initial condition
    P0 = 0; % Initial pressure, e.g., in Pa

    opts = odeset('RelTol',1e-6,'AbsTol',1e-8);

    tspan = [min(FlowTime) max(FlowTime)];

    % Solve ODE
    [tSol, PSol] = ode45(odefun, tspan, P0, opts);

end

%% WK3 Function
function [tSol, PSol] = WK3(Rp, Rd, C, Pd, FlowRate, FlowTime)

    % Interpolation of Q(t)
    Q_fun = @(t) interp1(FlowTime, FlowRate, t, 'linear', 'extrap');
    
    % Numerical differentiation of Q(t)
    dQdt = gradient(FlowRate, FlowTime);

    dQdt_fun = @(t) interp1(FlowTime, dQdt, t);

    % Define the differential equation
    % C * dP/dt = (1 + Rp/Rd) * Q(t) + Rp * C * dQ(t)/dt - (P - Pd) / Rd
    odefun = @(t, P) ((1 + Rp / Rd) * Q_fun(t) + Rp * C * dQdt_fun(t) - (P - Pd) / Rd) / C;

    % Initial condition
    P0 = 0; % Initial pressure

    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

    tspan = [min(FlowTime) max(FlowTime)];

    % Solve ODE
    [tSol, PSol] = ode45(odefun, tspan, P0, opts);

end

%% WK4Series Function
function [tSol, PSol] = WK4Series(Rp, Rd, C, L, Pd, FlowRate, FlowTime)

    % Interpolation of Q(t)
    Q_fun = @(t) interp1(FlowTime, FlowRate, t, 'linear', 'extrap');
    
    % First derivative of Q(t)
    dQdt = gradient(FlowRate, FlowTime);
    dQdt_fun = @(t) interp1(FlowTime, dQdt, t, 'linear', 'extrap');
    
    % Second derivative of Q(t)
    d2Qdt2 = gradient(dQdt, FlowTime);
    d2Qdt2_fun = @(t) interp1(FlowTime, d2Qdt2, t, 'linear', 'extrap');

    % Define the differential equation
    % C * dP/dt = (1 + Rp/Rd) * Q(t) + (L/Rd + Rp * C) * dQ(t)/dt + C * L * d2Q(t)/dt^2 - (P - Pd) / Rd
    odefun = @(t, P) ...
        ((1 + Rp / Rd) * Q_fun(t) + ...
         (L / Rd + Rp * C) * dQdt_fun(t) + ...
         C * L * d2Qdt2_fun(t) - ...
         (P - Pd) / Rd) / C;

    % Initial condition
    P0 = 0; % Initial pressure

    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

    tspan = [min(FlowTime) max(FlowTime)];

    % Solve ODE
    [tSol, PSol] = ode45(odefun, tspan, P0, opts);

end

%% WK4Parallel Function
function [tSol, PSol] = WK4Parallel(Rp, Rd, C, L, Pd, FlowRate, FlowTime)

    % Interpolation of Q(t)
    Q_fun = @(t) interp1(FlowTime, FlowRate, t, 'linear', 'extrap');
    
    % First derivative of Q(t)
    dQdt = gradient(FlowRate, FlowTime);
    dQdt_fun = @(t) interp1(FlowTime, dQdt, t, 'linear', 'extrap');
    
    % Second derivative of Q(t)
    d2Qdt2 = gradient(dQdt, FlowTime);
    d2Qdt2_fun = @(t) interp1(FlowTime, d2Qdt2, t, 'linear', 'extrap');
    
    % Define the differential equation (converted to a system of first-order ODEs)
    odefun = @(t, P) [
        P(2); % dP1/dt = P2 (pressure rate of change)
        P(3); % dP2/dt = P3 (second derivative of pressure)
        % Use the equation: (C*L/Rp)*d2P/dt^2 = Q(t) + L*((Rp + Rd)/(Rp*Rd))*dQ/dt + C*L*d2Q/dt^2 - (C + L/(Rp*Rd))*dP/dt - (P - Pd)/Rd
        (Q_fun(t) + L * ((Rp + Rd) / (Rp * Rd)) * dQdt_fun(t) + C * L * d2Qdt2_fun(t) - ...
         (C + L / (Rp * Rd)) * P(2) - (P(1) - Pd) / Rd) * (Rp / (C * L));
    ];
    
    % Initial conditions for the system
    P0 = 0;       % Initial pressure P(t)
    dPdt0 = 0;    % Initial first derivative of pressure
    d2Pdt20 = 0;  % Initial second derivative of pressure
    
    % Set options for ODE solver
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

    tspan = [min(FlowTime) max(FlowTime)];

    % Solve the system of ODEs using ode45
    [tSol, PSol] = ode15s(odefun, tspan, [P0, dPdt0, d2Pdt20], opts);

end

%% Coronary LPN Function
function [tSol, PSol] = CoronaryLPN(Ra, Ram, Rv, Rvm, Ca, Cim, Pv, PimScaling, Pim, PimTime, FlowRate, FlowTime)

    % Interpolation of Q(t)
    Q_fun = @(t) interp1(FlowTime, FlowRate, t, 'linear', 'extrap');
    
    % First derivative of Q(t)
    dQdt = gradient(FlowRate, FlowTime);
    dQdt_fun = @(t) interp1(FlowTime, dQdt, t, 'linear', 'extrap');

    % Interpolation of Pim(t)
    ScaledPim = PimScaling * Pim;

    Pim_fun = @(t) interp1(PimTime, ScaledPim, t, 'linear', 'extrap');
    
    % First derivative of Pim(t)
    dPimdt = gradient(ScaledPim, PimTime);
    dPimdt_fun = @(t) interp1(PimTime, dPimdt, t, 'linear', 'extrap');

    % Define the system of 2 ODEs
    odefun = @(t, P) [
        
        (1 / Ca) * (Q_fun(t) * (1 + Ra / Ram) + Ca * Ra * dQdt_fun(t) - (P(1) - P(2)) / Ram);
    
        (1 / Cim) * (Q_fun(t) - Ca * gradient(P(1),t) - (P(2) - Pv) / (Rv + Rvm) + Cim * dPimdt_fun(t))
    ];

    % Initial conditions
    P0 = [0, 0];
    
    % Set options for ODE solver
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

    tspan = [min(FlowTime) max(FlowTime)];

    % Solve the system of ODEs using ode15s
    [tSol, PSol] = ode15s(odefun, tspan, P0, opts);

end