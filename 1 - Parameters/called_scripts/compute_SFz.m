function SF = compute_SFz(t, Vm, Iion, Im)

nCells = size(Vm, 1);
if length(t)/t(end) > 5e5
    next_shift = 1e4;
else
    next_shift = 100;
end

APs = findpeaks(Vm(1, :), 'MinPeakHeight', 0.5 * max(Vm(1, :)), 'MinPeakDistance', next_shift);
Qion = zeros(nCells, length(APs));
Qm = zeros(nCells, length(APs));

for i = 1:nCells % cycle frontier cells
    
    [Iion_max, Iion_max_pos] = findpeaks(abs(Iion(i, :)), 'MinPeakHeight', 0.5 * max(abs(Iion(i,:))), 'MinPeakDistance', next_shift);
    [Im_max, Im_max_pos] = findpeaks(abs(Im(i, :)), 'MinPeakHeight', 0.5 * max(abs(Im(i,:))), 'MinPeakDistance', next_shift);
    idx_start_A = [Iion_max_pos - next_shift/2, 0];
    idx_start_B = [Im_max_pos - next_shift/2, 0];

    if ~isempty(Iion_max_pos)
        for j = 1:length(Im_max_pos) % cycle APs
            
            t1A = find(Iion(i, idx_start_A(j):end) < -0.01*Iion_max(j), 1, 'first') + idx_start_A(j)-2;
            t1B = find(Im(i, idx_start_B(j):end) < -0.01*Im_max(j), 1, 'first') + idx_start_B(j)-2;
            
            t2A = find( (t(t1A:end-1) > t(Iion_max_pos(j)))' & (Iion(i, t1A:end-1) > 0), 1, 'first') + t1A -1;
            t2B = find( (t(t1B:end-1) > t(Im_max_pos(j)))' & (Im(i, t1B:end-1) > 0), 1, 'first') + t1B -1;
            
            %%{
            subplot(311)
            plot(t, Vm)
            subplot(312)
            plot(t, Iion)
            hold on
            plot(t(t1A), Iion(i, t1A), 'r*')
            plot(t(t2A), Iion(i, t2A), 'go')
            subplot(313)
            plot(t, Im)
            hold on
            plot(t(t1B), Im(i, t1B), 'r*')
            plot(t(t2B), Im(i, t2B), 'go')            
            
            %}         
            
            Qion(i,j) = trapz(t(t1A:t2A), Iion(i, t1A:t2A));
            Qm(i,j) = trapz(t(t1B:t2B), Iion(i, t1B:t2B));   
            
            idx_start_A(j+1) = t2A + next_shift;
            idx_start_B(j+1) = t2B + next_shift;

        end
    end
    
    clc
    disp(['Computing SF for Cell #', num2str(i)])

end

SF = mean(Qion./Qm, 2);

end
