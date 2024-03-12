close all
clear
clc

model_name = 'koivu';

%% Settings
if strcmpi(model_name, 'TTP')
    scale_ms = 1e3;
    Cm = 0.185; % uF
    fcn = @tentusscher_panfilov_2006_a;
    Y0 = [2.72338164826056e-05;0.999590729087436;0.999992451362514;0.999372781405661;0.406199160609721;2.54387797721408e-05;0.000126461287514218;0.994539177919273;0.786133270192846;0.786127132867030;0.00121921087122944;-86.8326800357191;142.882218422513;0.000168170278947230;0.487842794634044;0.00288588409116802;3.35502269029023;1.85009268247943e-08;0.999998433285328];
    V_idx = 12;
    iStim = -(0.5:0.5:50); % pA
elseif strcmpi(model_name, 'koivu')
    scale_ms = 1;
    Cm = 0.05; %nF
    fcn = @koivumaki_2011;
    load('koivu_500s_ctrl.mat')
    Y0 = Y_end;
    V_idx = 24;
    iStim = -[0:0.01:0.1, (0.1:0.1:10)]*1e3; % pA
else
    disp('Error: choose appropriate model')
end

%% Stimulus
tStim = [0.1, 0.2, 0.3, 0.4, 1, 2, 3, 4, 10, 20, 30, 40, 100, 200, 300, 400]*1e-3*scale_ms; % ms , 10, 20, 30, 40, 100

qStim = tStim'.*iStim; % fC
Vthr  = -48; %-48;

%% Model
% Y0 = [0.59984, 0.64913, 0.21603, 0.00205, 0.68492, 0.41837, 6.0e-5, 0.5753, 0.39871, 0.57363, 3.0e-5, 0.99981, 4.6e-4, 0.30752, 1.6e-4, 0.76898, 0.02032, 0.02981, 0.01442, 0.23532, 0.67476, 7.305e-5, -69.83663, 0.706, 0.61493, 0.01309];

t_start = 0;
t_end = 1*scale_ms;

options = odeset('MaxStep', 1e-4*scale_ms, 'RelTol', 1e-3, 'AbsTol', 1e-6, 'OutputFcn',@myoutput); %% odewbar

%% Simulate
beat_mat = zeros(size(qStim));
time = cell(size(qStim));
Vm = cell(size(qStim));
Qall = zeros(length(tStim), length(iStim));

for j = 1:length(tStim)
    for k = 1:length(iStim)
        [time{j, k}, Y] = ode15s(@(x,y) fcn(x, y, tStim(j), iStim(k)), [t_start, t_end], Y0, options);
        %       [time{j, k}, Y] = ode15s(@(x,y) lindblad_fixed_Nai_1996(x, y, 1, tStim(j), iStim(k)), [t_start, t_end], Y0, options);
        
        Vm{j, k} = Y(:, V_idx);
        
        clc
        disp(['Simulating time ... ', num2str(j), '/', num2str(length(tStim)), ' and stimulus...', num2str(k), '/', num2str(length(iStim))])
        
        if isreal(Vm{j,k}) && ~isempty((Vm{j,k}))
            [OS, OSpos] = findpeaks(Vm{j, k}, 'MinPeakHeight', 0);
            if isempty(OSpos)
                beat_mat(j, k) = 0;
                Qall(j, k) = nan;
            else
                beat_mat(j, k) = 1;
                
                %%% Qthr
                Inet = -diff(Vm{j,k}*1e-3) ./ diff(time{j,k}); % = [F*V/s]= [A]                    
                %%% 1: Trovo i picchi della corrente (gli istanti)
                [Imax, t_Imax] = findpeaks(-Inet, 'MinPeakHeight', max(-Inet)/5);
                              
                i_stimulus = zeros(1, length(time{j,k}));
                i_stimulus(time{j,k} < tStim(j)) = iStim(k);
                
                %%% 2: Prendo l'istante a 1% del picco -> t1%
                if ~isempty(t_Imax)
                    t1 = 1; %find(abs(Inet) > 0.01*Imax(1), 1, 'first');
                    %%% 3: Istante in cui Im = 0
                    %                 t2 = find((time{j, k}(t1:end-1) > time{j,k}(t_Imax(1))) & (Inet(t1:end) > 0), 1, 'first') + t1;
                    t2 = find( Inet(t1:end) > 0, 1, 'first') + t1;
                    
                    %             t2 = find(Vm{j,k}(t1:end-1) > Vthr, 1, 'first')+t1; % time to Vthr = -48 mV
                    
                    Inet = i_stimulus(1:end-1);
                    Inet_int = Inet(t1:t2);
                    %             Inet_int(Inet_int > 0) = 0;
                    
                    Qall(j, k) = trapz(time{j,k}(t1:t2), Inet_int); % pC/pF

                else
                    Qall(j, k) = nan;
                end
                if exist('t1', 'var')
                    plot(time{j,k}(1:end-1)/scale_ms, Inet)
                    hold on
                    plot(time{j,k}(t1)/scale_ms, Inet(t1), 'r*')
                    plot(time{j,k}(t2)/scale_ms, Inet(t2), 'go')
                    xlabel('t [s]'), ylabel('nA')
                end
            end
            
        else
            Qall(j, k) = nan;
        end
    end
end

%% Plot Qhtr
tStim = tStim*1e3/scale_ms;

figure
for i = 1:size(Qall, 2)
    plot(tStim, Qall(:, i), 'k*')
    hold on
    
    if i == 1
        xlim([0 max(tStim)*1.2])
%         ylim([-1 0])
    end
end
xticks(tStim)
xticklabels(tStim)
xlabel('Stimulus time [ms]')
ylabel('Q_{thr} [pC/pF]')

Qthr = -min(abs(Qall), [], 2, 'omitnan');

plot(tStim , Qthr, 'ro')

%% Qall map
figure
imagesc(Qall, 'AlphaData', ~isnan(Qall))
xticks(1:length(iStim))
yticks(1:length(tStim))
xticklabels(-iStim/1e3)
yticklabels(tStim)
xlabel('I_{stim} [nA]')
ylabel('t_{stim} [ms]')
title('Q_{thresh}')

%% Linear regression
Qthr_2pA = Qthr; %(:,4);
coeff = polyfit(tStim(~isnan(Qthr_2pA)), Qthr_2pA(~isnan(Qthr_2pA))', 1);
yy = tStim*coeff(1) + coeff(2);

figure
plot(tStim, Qthr_2pA, 'k*')
hold on
plot(tStim, yy)
xlabel('Time [ms]'), ylabel('Q_{thr} [pC/pF]')

%% Conversion to nC/cm^2 (assuming specific capacitance = 1 uF/cm^2)
% Qthr_on_cm2 = Qthr * 1000; % [nC/nF]*[1e-3 uF/cm^2] = [nC/cm^2] %%% for Inet
Qthr_on_cm2 = Qthr * (1/Cm); % [pC]*[1/pF]*[1 uF/cm^2] = [pC/pF]* [1e6 pF/cm^2]= 1e3 [nC/cm^2] %%% for iStim

disp('            [nC/cm^2]')
disp(Qthr_on_cm2')

%% Save
save(['Qall_', model_name, '.mat'], 'Qall')
save(['stimResults_Vm_', model_name, '.mat'], 'time', 'Vm')
save(['coeff_Qthr_', model_name,'.mat'], 'coeff')

%%

function stop = myoutput( x, optimValues, state ) %#ok<INUSL>
%MYOUTPUT stop intergration based on some criterion
stop = false;
persistent numMajorTimeStep
global numDiffEqEval
switch state
    case []
        % several output values may be calculated in one major time
        % step
        numMajorTimeStep = numMajorTimeStep + numel(x);
        
        if numMajorTimeStep > 20000 || numDiffEqEval > 5000
            stop = true;
        end
    case 'init'
        % Setup for plots or guis
        numMajorTimeStep = 0; numDiffEqEval = 0;
    case 'done'
        % Cleanup of plots, guis, or final plot
end
end
