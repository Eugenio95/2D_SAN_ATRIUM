function CV = CV_in_SEP(t, Vm_mat, x_s1, y_s1, cell_length, nCells)

ind_start = sub2ind([200, 200], x_s1, y_s1);
ind_end   = sub2ind([200, 200], x_s1, y_s1 + nCells);

% figure
% plot(t, Vm_mat(ind_start, :))
% hold on
% plot(t, Vm_mat(ind_end, :))

[~, OSpos_start] = findpeaks(diff(Vm_mat(ind_start, :))./diff(t)', 'MinPeakHeight', 3e3, 'MinPeakDistance', 50);
[~, OSpos_end] = findpeaks(diff(Vm_mat(ind_end, :))./diff(t)', 'MinPeakHeight', 3e3, 'MinPeakDistance', 50);

if length(OSpos_start) > 1
    if OSpos_start(1) < OSpos_end(1)
        n_peak_min = min(length(OSpos_start), length(OSpos_end));
        delta_t = t(OSpos_end(1:n_peak_min)) - t(OSpos_start(1:n_peak_min));
    else
        n_peak_min = min(length(OSpos_start), length(OSpos_end))-1;
        delta_t = t(OSpos_end(2:n_peak_min)) - t(OSpos_start(1:n_peak_min-1));
    end
    
    delta_s = nCells * cell_length;
    CV = (delta_s ./ delta_t) * 100; % (# of cell * length of cells)/dt [cm/s]
    disp(['CV = ', num2str(round(mean(CV), 2)), '+', num2str(round(std(CV), 2)), ' cm/s', ' || t_SEP = ', num2str(round(mean(delta_t), 2)), '+', num2str(round(std(delta_t), 2)), ' s'])
    
else
    CV = nan;
    disp('Cells do not beat')
end

end

