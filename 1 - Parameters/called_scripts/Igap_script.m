
idx_atrio = 0;
idx_san = 1;
idx_fibro = 3;
idx_fat = 9;

% Creo file capacità
Cm = geom_mat;
Cm(Cm == idx_atrio) = 50e-12;
Cm(Cm == idx_san) = 32e-12;
Cm(Cm == idx_fibro) = 12.4e-12;
Cm(Cm == idx_fat) = inf;

%% Calcolo Igap
Igap   = zeros(200*200, length(t));

disp('Inizio calcolo Igap')
for i = 1:size(geom_mat, 1)
    for j = 1:size(geom_mat, 2)
               
        ind_ij = sub2ind(size(Cm), i, j);
        
        if j > 1
            ind_j_minus_1 = sub2ind(size(Cm), i, j-1);
            Igap_L = squeeze((Vm_mat(ind_j_minus_1, :) - Vm_mat(ind_ij, :)) * gJ(i, j, 1))';
        else
            Igap_L = squeeze(zeros(1, size(Vm_mat, 3)));
        end
        
        if i > 1
            ind_i_minus_1 = sub2ind(size(Cm), i-1, j);
            Igap_U = squeeze((Vm_mat(ind_i_minus_1, :) - Vm_mat(ind_ij, :)) * gJ(i, j, 2))';
        else
            Igap_U = squeeze(zeros(1, size(Vm_mat, 3)));
        end
        
        if j < size(geom_mat, 2)
            ind_j_plus_1 = sub2ind(size(Cm), i, j+1);
            Igap_R = squeeze((Vm_mat(ind_j_plus_1, :) - Vm_mat(ind_ij, :)) * gJ(i, j, 3))';
        else
            Igap_R = squeeze(zeros(1, size(Vm_mat, 3)));
        end
        
        if i < size(geom_mat, 1)
            ind_i_plus_1 = sub2ind(size(Cm), i+1, j);
            Igap_D = squeeze((Vm_mat(ind_i_plus_1, :) - Vm_mat(ind_ij, :)) * gJ(i, j, 4))';
        else
            Igap_D = squeeze(zeros(1, size(Vm_mat, 3)));
        end
    
        % Negative sign: incoming current is depolarizing
        Igap(ind_ij,:) = -(Igap_L + Igap_U + Igap_R + Igap_D) * 1e-3/ Cm(i, j); % in A/F
    
    end        
    
    clc
    disp(['Igap column = ', num2str(i), '/200'])
end
clear Igap_L Igap_U Igap_R Igap_D
% Igap = Igap * 1e-3 *50; % per qualche cazzo di motivo ci vuole, sennò non è paragonabile a Inet

clc
disp('Fine calcolo Igap')
% Calcolo corrente totale cellule da variazione Vm (devo togliere Igap)
% "diff" fa la differenza in avanti
Inet = -(diff(Vm_mat*1e-3, 1, 2))./diff(t)'; % A/F      %%% A/F = C/V * V/s * (1/(C/V)) 

%% Plot
figure
subplot(221)
plot(t(1:end-1), squeeze(Inet(1900,:)))
title('I_{tot} atrium')
ylabel('pA/pF')
subplot(222)
plot(t, squeeze(Igap(1900,:)))
title('I_{gap} atrium')

subplot(223)
plot(t(1:end-1), squeeze(Inet(91*200+99,:)))
title('I_{tot} SAN')
ylabel('pA/pF')
subplot(224)
plot(t, squeeze(Igap(91*200+ 99,:)))
title('I_{gap} SAN')
