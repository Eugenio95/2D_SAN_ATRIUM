function SF = compute_SF(t, Vm_mat, Inet, Iion, Igap, nCells, flagStim, iStim_amp)

thr = load('coeff_Qthr.mat');

xx = linspace(0, 4e-3, 1e3);
yy = thr.coeff(1)*xx + thr.coeff(2);

APs = findpeaks(abs(Inet(1, :)), 'MinPeakHeight', 0.5 * max(abs(Inet(1, :))), 'MinPeakDistance', 100);

iStim = zeros(1, length(Inet));
for k = 0:numel(APs)-1
    iStim( find(t > 0.01+k*(1/3), 1, 'first') : find(t > 0.01+k*(1/3)+2e-3, 1, 'first')) = iStim_amp/50e-3; % pA/pF
end

if length(t)/t(end) > 5e5
    next_shift = 1e4;
else
    next_shift = 100;
end

% figure
% plot(xx, yy)
% xlabel('t [s]'), ylabel('pC')

Qgap = zeros(nCells, numel(APs));
Qthr = zeros(nCells, numel(APs));
Qstim = zeros(nCells, numel(APs));
for i = 1:nCells % cycle frontier cells
    
    [Inet_max, Inet_max_idx] = findpeaks( -Inet(i, :), 'MinPeakHeight', 0.5*max(-Inet(i, :)), 'MinPeakDistance', next_shift); % on Itot
%     [~, Itot_max_idx] = findpeaks( -Inet(i, :), 'MinPeakHeight', 0.5*max(-Inet(i, :)), 'MinPeakDistance', 5000); % on Itot
    idx_start = [Inet_max_idx- 30, 0];
    
    if ~isempty(Inet_max_idx)
        for j = 1:length(Inet_max_idx) % cycle APs
            
            t1 = find(Inet(i, idx_start(j):end) < -0.01*Inet_max(j), 1, 'first') + idx_start(j)-2;
            %%% 3: Istante in cui Im = 0
            t2 = find(t(t1:end-1) > t(Inet_max_idx(j)) & (Inet(i, t1:end-1) > 0), 1, 'first') + t1 -1;
            
            Qgap(i,j)  = trapz(t(t1:t2), Igap(i, t1:t2)); % if I in pA/pF => Q in pC/pF
            tA = t(t2)-t(t1);
                        
            %%% Trucco per evitare intervalli troppo lunghi quando non ho
            %%% ho il notch... legittimo...??
%             if tA > 0.004
%                 tA = 0.004;
%             end    
                        
            Q_estim_idx = find(abs(xx - tA) < 1e-4, 1, 'first');
            if isempty(Q_estim_idx)
                Q_estim_idx = find(abs(xx - tA) < 1e-3, 1, 'first');
                if isempty(Q_estim_idx)
                    Q_estim_idx = find(abs(xx - tA) < 1e-2, 1, 'first');
                    if isempty(Q_estim_idx)
                        Qthr(i, j) = nan;
                    else
                        Qthr(i, j) = yy( Q_estim_idx ); % pC
                    end
                else
                    Qthr(i, j) = yy( Q_estim_idx ); % pC
                end
            else
                Qthr(i, j) = yy( Q_estim_idx ); % pC
            end
            
            if flagStim(i) == 1
                t2_stim = find(t > t(t1) + 2e-3, 1, 'first'); % only stimulus duration
                Qstim(i, j) = trapz(t(t1:t2_stim), iStim(t1:t2_stim)); % pA/pF * s = pC/pF
            else
                Qstim(i, j) = 0;
            end
            
            idx_start(j+1) = t2 + next_shift;
       
%{             
            sgtitle(['Cell ', num2str(i)])
            subplot(511)
            plot(t, Vm_mat(i, :))
            
            subplot(512)
            plot(t, Inet(i, :))
            hold on
            plot(t(t2), Inet(i, t2), 'ko')
            plot(t(t1), Inet(i, t1), 'r*')
            ylabel('pA/pF')
            title('Inet')
            hold off
            
            subplot(513)
            plot(t, Iion(i, :))
            hold on
            plot(t(t2), Iion(i, t2), 'ko')
            plot(t(t1), Iion(i, t1), 'r*')
            ylabel('pA/pF')
            title('Itot')
            hold off
            
            subplot(514)
            plot(t, Igap(i, :))
            hold on
            plot(t(t2), Igap(i, t2), 'ko')
            plot(t(t1), Igap(i, t1), 'r*')
            ylabel('pA/pF')
            title('Igap')
            hold off
            
            subplot(515)
            plot(t, iStim)
            hold on
            plot(t(t1), iStim(t1), 'r*')
            plot(t(t2), iStim(t2), 'ko')
            ylabel('pA/pF')
            title('iStim')
            hold off
%}
            
        end
    else
        for j = 1:length(Inet_max_idx) % cycle APs
            Qgap(i,j) = nan;
        end
    end
    
    clc
    disp(['Computing Qgap for cells at frontier...', num2str(round(i/nCells * 100)), '%'])
    
end

Qthr = -0.02; % it is almost constant vs time if we use Vthr

% SF = round(mean( ((Qgap) - (Qstim)) ./ (Qthr), 2, 'omitnan'), 1);
SF = mean( ((Qgap) - (Qstim)) ./ (Qthr), 2, 'omitnan');

end

