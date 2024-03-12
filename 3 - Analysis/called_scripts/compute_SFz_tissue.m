function SF = compute_SFz_tissue(cell_idx, t, Vm_mat, Iion, Im, geom_mat)

if length(t)/t(end) > 5e5
    next_shift = 1e4;
else
    next_shift  = 2000; % about half the number of samples between two AP
    first_shift = 1000;
end

%% Calcolo carica di gap junctions
APs = findpeaks(Vm_mat(cell_idx(1), :), 'MinPeakHeight', 0.5 * max(Vm_mat(cell_idx(1), :)), 'MinPeakDistance', next_shift);
Qion = zeros(length(cell_idx), length(APs));
Qm = zeros(size(Qion));

for i = 1:length(cell_idx) % cycle frontier cells
    
    if geom_mat(cell_idx(i)) == 0 || geom_mat(cell_idx(i)) == 1
        
        [Iion_max, Iion_max_pos] = findpeaks( -(Iion(cell_idx(i), :)), 'MinPeakHeight', 0.1*max(-(Iion(cell_idx(i), :))), 'MinPeakDistance', next_shift);
        [Im_max, Im_max_pos] = findpeaks( -(Im(cell_idx(i), :)), 'MinPeakHeight', 0.1 * max(-(Im(cell_idx(i),:))), 'MinPeakDistance', next_shift);
        
        Iion_max_pos(Iion_max_pos < first_shift | Iion_max_pos > (length(t) - next_shift)) = [];
        Im_max_pos(Im_max_pos < first_shift | Im_max_pos > (length(t) - next_shift)) = [];
        
        idx_start_A = [Iion_max_pos - next_shift+1, 0];
        idx_start_B = [Im_max_pos - next_shift+1, 0];
        
        if ~isempty(Iion_max_pos)
            for j = 1:length(Iion_max_pos) % cycle APs
                
                try
                    t1A = find((Iion(cell_idx(i), idx_start_A(j):end) < -0.01*Iion_max(j)) & (Iion(cell_idx(i), idx_start_A(j):end) < 0), 1, 'first') + idx_start_A(j);
                    t1B = find((Im(cell_idx(i), idx_start_B(j):end) < -0.01*Im_max(j)) & (Im(cell_idx(i), idx_start_B(j):end) < 0), 1, 'first') + idx_start_B(j);
                    
                    t2A = find( (t(t1A:end-1) > t(Iion_max_pos(j)))' & (Iion(cell_idx(i), t1A:end) > 0), 1, 'first') + t1A;
                    t2B = find( (t(t1B:end-1) > t(Im_max_pos(j)))' & (Im(cell_idx(i), t1B:end) > 0), 1, 'first') + t1B;
                    
                    Qion(i,j) = trapz(t(t1A:t2A), Iion(cell_idx(i), t1A:t2A));
                    Qm(i,j) = trapz(t(t1B:t2B), Im(cell_idx(i), t1B:t2B));
%                     Qm(i,j) = trapz(t(t1A:t2A), Im(cell_idx(i), t1A:t2A));
                    
                    idx_start_A(j+1) = t2A + next_shift;
                    idx_start_B(j+1) = t2B + next_shift;
                    
                catch
                    Qion(i,j) = nan;
                    Qm(i,j)   = nan;
                end
                
            end
            
        else
            for j = 1:length(Iion_max_pos) % cycle APs
                Qion(i,j) = nan;
                Qm(i,j)   = nan;
            end
        end
        
    else
        for j = 1:length(Iion_max_pos) % cycle APs
            Qion(i,j) = nan;
            Qm(i,j)   = nan;
        end
    end
    
    if mod(i, 100) == 0
        
        %{
            subplot(311)
            plot(t, Vm_mat(cell_idx(i), :))
            hold on
            plot(t(t1A), Vm_mat(cell_idx(i), t1A), 'r*')
            plot(t(t2A), Vm_mat(cell_idx(i), t2A), 'go')            
            
            subplot(312)
            plot(t(1:end-1), Iion(cell_idx(i), :))
            hold on
            plot(t(t1A), Iion(cell_idx(i), t1A), 'r*')
            plot(t(t2A), Iion(cell_idx(i), t2A), 'go')
            title('Iion')
            hold off
            
            subplot(313)
            plot(t(1:end-1), Im(cell_idx(i), :))
            hold on
            plot(t(t1B), Im(cell_idx(i), t1B), 'r*')
            plot(t(t2B), Im(cell_idx(i), t2B), 'go')
            hold off
            title('Im')
        %}
        
        clc
        disp(['Computing Qgap for cells at frontier...', num2str(round(i/length(cell_idx) * 100)), '%'])
    end
    
end

%% Calcolo Safety Factor
SF = Qion./Qm;

end
