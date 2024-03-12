function TOPpos = find_TOP_Kohadja2020(Vm)

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%% SINGLE CELL %%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% TOP = 15% dV^2/dt^2_max
%%% from Kohadja et al. 2020 

Vm = squeeze(Vm);

% APD -> 15% (dV^2/dt^2)_max
dV2_dt2_max = diff(diff(squeeze(Vm)));
[TOP, pos_temp] = findpeaks(dV2_dt2_max, 'MinPeakHeight', 0.25*max(dV2_dt2_max), 'MinPeakDistance',100);
% con findpeaks trovo i massimi della derivata seconda e le
% loro posizioni

if length(TOP) > 2
    
    for i_top = 1:length(TOP)
        % prendo i punti in cui la derivata seconda è sopra soglia
        pos_above_thresh = find( diff(diff(squeeze(Vm))) > 0.15*TOP(i_top))';
        % prendo solo il primo di ogni gruppetto (uno per ogni battito)
        out_pos = pos_above_thresh(diff([0 pos_above_thresh]) > 10)';
        
        pos15_idx = 0;
        % controllo che gli istanti del TOP trovati siano
        % vicini al picco max della derivata seconda, nel caso
        % li salvo in un nuovo vettore (pos15_idx)
        for idx_out = 1:length(out_pos)
            
            if abs(out_pos(idx_out) - pos_temp(i_top)) < 50
                pos15_idx = [pos15_idx, out_pos(idx_out)];
            end
            
        end
        % tengo solo l'istante giusto (tolgo lo zero preallocato
        % ed eventuali istanti aggiuntivi che non sono stati scartati)
        if length(pos15_idx) > 2
            pos15_long = pos15_idx(pos15_idx ~= 0);
            pos15(i_top) = pos15_long(1);
        elseif length(pos15_idx) == 2
            pos15(i_top) = pos15_idx(pos15_idx ~= 0);
        else
            pos15(i_top) = pos15_idx;
        end
        
    end
    % salvo l'istante di TOP
    TOPpos = pos15(pos15 ~= 0);
else
    TOPpos = [];    
end

end

