%% Estraggo gli indici di tutte le cellule alla frontiera
% figure
% imagesc( geom_mat - geom_mat(:,[2:end end]) )
front_SEP_cells = find( geom_mat - geom_mat(:,[2:end end])  == thresh_frontier) +200; % find the atrial cells of the frontier

%% Coefficienti Qthr
addpath 2_cellule
thr = load('coeff_Qthr.mat');

xx = linspace(0, 100e-3, 1e3);
yy = thr.coeff(1)*xx + thr.coeff(2);

next_shift = 200;

%% Calcolo carica di gap junctions
Qgap = zeros(length(front_SEP_cells), 5);
Qthr = size(Qgap);

% aggiungo cellule a caso nell'atrio per testare SF
front_SEP_cells = [front_SEP_cells', 201, 7940, 20000, 31840, 39840];

for i = 1:length(front_SEP_cells) % cycle frontier cells
        
    [Inet_max, Inet_max_idx] = findpeaks( -Inet(front_SEP_cells(i), :), 'MinPeakHeight', 0.2*max(-Inet(front_SEP_cells(i), :)), 'MinPeakDistance', next_shift); % on Itot
    Inet_max = -Inet_max;
    idx_start = [Inet_max_idx - 30, 0];
    
    if ~isempty(Inet_max_idx)
        for j = 1:length(Inet_max_idx) % cycle APs
            
            t1 = find(Inet(front_SEP_cells(i), idx_start(j):end) < 0.01*Inet_max(j), 1, 'first') + idx_start(j)-2;
            %%% 3: Istante in cui Im = 0
            t2 = find(t(t1:end-1)' > t(Inet_max_idx(j)) & Inet(front_SEP_cells(i), t1:end) > 0, 1, 'first') + t1-1;
            
            Igap_int = Igap(front_SEP_cells(i), t1:t2);
            Qgap(i,j)  = trapz(t(t1:t2), Igap_int); % pC
            tA = t(t2)-t(t1);
%             disp(['tA = ', num2str(tA), 's'])
            
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
            
            idx_start(j+1) = t2 + next_shift;
            
%%{
            subplot(311)
            plot(t, Vm_mat(front_SEP_cells(i), :))
            
            subplot(312)
            plot(t(1:end-1), Inet(front_SEP_cells(i), :))
            hold on
            plot(t(t2), Inet(front_SEP_cells(i), t2), 'ko')
            plot(t(t1), Inet(front_SEP_cells(i), t1), 'r*')
            title('Inet')
            hold off
            
            subplot(313)
            plot(t, Igap(front_SEP_cells(i), :))
            hold on
            plot(t(t2), Igap(front_SEP_cells(i), t2), 'ko')
            plot(t(t1), Igap(front_SEP_cells(i), t1), 'r*')
            title('Igap')
            hold off
            
            %}
        end
    else
        for j = 1:length(Inet_max_idx) % cycle APs
            Qgap(i,j) = nan;
        end
    end
    
    clc
    disp(['Computing Qgap for cells at frontier...', num2str(round(i/length(front_SEP_cells) * 100)), '%'])
    
end

Qthr = -0.02; % it is almost constant vs time if we use Vthr

%% Calcolo Safety Factor
SF = mean(Qgap./Qthr, 2, 'omitnan');

disp(['SF of central cell = ', num2str(SF(25)) ])
disp(['Mean SF = ', num2str(round(mean(SF))) ])
