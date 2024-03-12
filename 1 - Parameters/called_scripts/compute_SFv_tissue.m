function SF = compute_SFv_tissue(cell_idx, t, Vm_mat, Inet, Igap, geom_mat)

%% Coefficienti Qthr
addpath called_scripts
thr = load('coeff_Qthr_koivu.mat');

xx = linspace(0, 4, 1e3);
yy = thr.coeff(1)*xx + thr.coeff(2);

if length(t)/t(end) > 5e5
    next_shift = 1e4;
else
    next_shift  = 2200; % about half the number of samples between two AP
    first_shift = 1000;
end

%% Calcolo carica di gap junctions
APs = findpeaks(Vm_mat(cell_idx(1), :), 'MinPeakHeight', 0.5 * max(Vm_mat(cell_idx(1), :)), 'MinPeakDistance', next_shift);
Qgap = zeros(length(cell_idx), length(APs));
Qthr = zeros(size(Qgap));
noBeat = 0;

for i = 1:length(cell_idx) % cycle frontier cells
       
    if geom_mat(cell_idx(i)) == 0
        %%% atrium
        dVdt = diff(Vm_mat(cell_idx(i), :))./diff(t');
        [dVdt_max, ~] = findpeaks(dVdt, 'MinPeakHeight', max(dVdt)/2, 'MinPeakDistance', next_shift);
        
        [Inet_max, Inet_max_idx] = findpeaks( -Inet(cell_idx(i), :), 'MinPeakHeight', 0.2*max(-Inet(cell_idx(i), :)), 'MinPeakDistance', next_shift); % on Itot
        Inet_max = -Inet_max;
        Inet_max_idx(Inet_max_idx < first_shift | Inet_max_idx > (length(t) - next_shift)) = [];
        idx_start = [Inet_max_idx - first_shift+1, 0];
        
        [V_max, V_max_idx] = findpeaks( Vm_mat(cell_idx(i), :), 'MinPeakHeight', 0, 'MinPeakDistance', next_shift); % on Itot
        V_max = -V_max;
        V_max_idx(V_max_idx < first_shift | V_max_idx > (length(t) - next_shift)) = [];
        V_idx_start = [V_max_idx - first_shift+1, 0];
        
%         [~, Igap_max_pos] = findpeaks(-Igap(cell_idx(i), :), 'MinPeakHeight', max(-Igap(cell_idx(i), :))/2, 'MinPeakDistance', next_shift);
        
        if ~isempty(V_max_idx)
            for j = 1:length(V_max_idx) % cycle APs
                
                try
                    V_rest = min(Vm_mat(cell_idx(i), V_idx_start(j):V_max_idx(j)));
                    t1b = find(Vm_mat(cell_idx(i), V_idx_start(j):end) > (0.05*(V_max(j)-V_rest)+V_rest), 1, 'first') + V_idx_start(j);
                    t1 = find(Inet(cell_idx(i), idx_start(j):end) < 0.01*Inet_max(j), 1, 'first') + idx_start(j);
                    %%% 3: Istante in cui Im = 0
                    %                     t2 = find(Inet(cell_idx(i), t1:end) > 0, 1, 'first') + t1-1;
                    if dVdt_max(j) < 1e5 || (t(t1)-t(t1b)) > 10e-3 %
                        t1 = find(Vm_mat(cell_idx(i), V_idx_start(j):end) > (0.01*(V_max(j)-V_rest)+V_rest), 1, 'first') + V_idx_start(j);
                        Igap_shift = round(next_shift/10);
                    else
                        Igap_shift = 0;
                        %                         t1 = find(Vm_mat(cell_idx(i), V_idx_start(j):end) > (0.01*(V_max(j)-V_rest)+V_rest), 1, 'first') + V_idx_start(j);
                        %                         t1 = find(Inet(cell_idx(i), idx_start(j):end) < 0.01*Inet_max(j), 1, 'first') + idx_start(j);
                        next_shift = 3000;
                    end
                    
                    Igap_idx_start = t1 - Igap_shift; %Igap_max_pos
                    t2 = find(Igap(cell_idx(i), Igap_idx_start:end) > 0, 1, 'first') + Igap_idx_start-2;
                    
                    Qgap(i,j)  = trapz(t(t1:t2), Igap(cell_idx(i), t1:t2)); % pC
                    tA = (t(t2)-t(t1))*1e3;
                    
                    if tA < 10
%                         Q_estim_idx = find(abs(xx - tA) < 1e-4, 1, 'first');
%                         if isempty(Q_estim_idx)
%                             Q_estim_idx = find(abs(xx - tA) < 1e-3, 1, 'first');
%                             if isempty(Q_estim_idx)
                                Q_estim_idx = find(abs(xx - tA) < 1e-2, 1, 'first');
                                if isempty(Q_estim_idx)
                                    Qthr(i, j) = nan;
                                else
                                    Qthr(i, j) = yy( Q_estim_idx ); % pC
                                end
%                             else
%                                 Qthr(i, j) = yy( Q_estim_idx ); % pC
%                             end
%                         else
%                             Qthr(i, j) = yy( Q_estim_idx ); % pC
%                         end
                    else
                        Qthr(i, j) = yy(end);
                    end
                    
                    V_idx_start(j+1) = t2 + next_shift;
                    
                catch
                    Qthr(i, j) = nan;
                end
                
                %                 check_plot_t1t2
                                
            end
        else
            noBeat = noBeat + 1;
            for j = 1:length(V_max_idx) % cycle APs
                Qgap(i,j) = nan;
            end
        end
        
    elseif geom_mat(cell_idx(i)) == 1
        %%% SAN
        
        
    else
        for j = 1:length(V_max_idx) % cycle APs
            Qgap(i,j) = nan;
        end
    end
    
    if mod(i, 400) == 0
        clc
        disp(['Computing Qgap for cells at frontier...', num2str(round(i/length(cell_idx) * 100)), '%'])
    end
    
end

%% Calcolo Safety Factor
SF = Qgap./Qthr; %%% mean(, 2, 'omitnan')
disp(['Number of non beating cells = ', num2str(noBeat)])

end
