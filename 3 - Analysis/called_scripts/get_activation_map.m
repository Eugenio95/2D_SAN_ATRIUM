function act_time_mat = get_activation_map(t, Vm_mat, geom_mat, front_SEP_cells)

for i = size(Vm_mat, 1):-1:1
    
    [~, dVdtmax_pos{i}] = findpeaks(diff(Vm_mat(i, :))./diff(t)', 'MinPeakHeight', 1e3, 'MinPeakDistance', 2000);
    
    %     if isempty(dVdtmax_pos{i})
    %         [~, dVdtmax_pos{i}] = findpeaks(diff(Vm_mat(i, :))./diff(t)', 'MinPeakHeight', 1e2, 'MinPeakDistance', 2000);
    %     end
    
    [~, wrnID] = lastwarn;
    warning('off', wrnID)
    
    if mod(i, 400) == 0
        clc
        disp(['Extracting TOPs for activation time... ', num2str( round((size(Vm_mat, 1)-i)/size(Vm_mat, 1)*100)), '%'])
    end
    
end

dVdtmax_pos_frontier = dVdtmax_pos(front_SEP_cells);

try
    front_act_pos = cell2mat(cellfun(@(x) x(end-2), dVdtmax_pos_frontier, 'UniformOutput', false));
    [zero_act_time, ~] = min(front_act_pos);
    
    for i = size(Vm_mat, 1):-1:1
        
        [j1, j2] = ind2sub([200, 200], i);
        
        if geom_mat(j1, j2) == 0 || geom_mat(j1, j2) == 1
            %         atrial_act_pos = dVdtmax_pos{i}( find(dVdtmax_pos{i} >= zero_act_time, 1, 'first') );
            %         act_time(i) = t(atrial_act_pos);
            %     elseif geom_mat(j1, j2) == 1
            %         try
            %             atrial_act_pos = dVdtmax_pos{i}( find(dVdtmax_pos{i} >= (zero_act_time - round(CLsan_min) +50), 1, 'first') );
            %             act_time(i) = t(atrial_act_pos);
            %         catch
            %             act_time(i) = nan;
            %         end
            try
                atrial_act_pos = dVdtmax_pos{i}( (find( abs(dVdtmax_pos{i} - zero_act_time) < min(diff(dVdtmax_pos{i}))/2, 1, 'first') ));
                act_time(i) = t(atrial_act_pos);
            catch
                act_time(i) = nan;
            end
        else
            act_time(i) = nan;
        end
        
        if mod(i, 400) == 0
            clc
            disp(['Computing activation times... ', num2str( round((size(Vm_mat, 1)-i)/size(Vm_mat, 1)*100)), '%'])
        end
    end
    act_time_mat = (reshape(act_time, 200, 200) - t(zero_act_time)) * 1e3; % in ms
    clc
    disp('Activation time computation concluded')

catch
    
    act_time_mat = nan(200); % in ms
    disp('Activation time computation concluded')
    
end


end
