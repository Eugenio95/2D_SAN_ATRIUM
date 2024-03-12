function g_rand = r_g_SEV(sigma, num_cells, c_flag, c_type, alpha, beta, atrial_tissue)

global x_c sep_high sep_half_width sep_low % y_c x_r y_r
    
% Geometry
Cm = 32;
L = 70;
R = 4;

% Nominal values
P_CaL = 0.2;   % nanoA_per_millimolar (in i_CaL)
P_CaT = 0.02;   % nanoA_per_millimolar (in i_CaT)
g_Kr = 0.0021637;   % microS (in i_Kr)
K_NaCa = 4.0;   % nanoA (in i_NaCa)
i_NaK_max = 0.063;   % nanoA (in i_NaK)
g_Na = 0.0125;   % microS (in i_Na)
g_to = 0.002;   % microS (in i_to)
g_Ks = 0.0016576; % microS (in i_Ks)
g_f = 0.03; % microS
g_KACh = 0.00864;   % microS (in i_KACh)
P_up_basal = 12.0;  % millimolar_per_second (in Ca_intracellular_fluxes)

%% Nominal (no )
% nomCond = [P_CaL, P_CaT, g_Kr, K_NaCa, i_NaK_max, g_Na, g_Ks, g_f, g_to, P_up_basal, g_KACh, Cm, L, R]; 
nomCond = [P_CaL, P_CaT, g_Kr, K_NaCa, i_NaK_max, g_Na, g_Ks, g_f, g_KACh, Cm, L, R]; 

% Per_factor = [11.4, 3.2, 5, 10, 6.5, 20, 5, 12.6, 1, 1, 1, 2.03, 1.23, 1.5]; %% From Ly & Weinberg (2018)
Per_factor = [11.4, 3.2, 5, 10, 6.5, 20, 5, 12.6, 1, 1, 1, 2.03, 1.23, 1.5]; %% From Ly & Weinberg (2018)
% where = 1 no data about that, but ony random hetero

%% Randomization
if c_flag == 0
    
    % Random hetero, no influence on space
    for i = num_cells:-1:1
        g_rand(i,:) = nomCond .* exp( sigma * randn(1,length(nomCond)) );
    end
    
elseif c_flag == 1
    
    if c_type == 0
        % Linear gradient from center to periphery
        for i = num_cells:-1:1
            g_scale_per = Per_factor/num_cells*i + 1 * (num_cells-i)/num_cells ;
            g_rand(i,:) = nomCond .* exp(sigma*randn(1,length(nomCond))) .* g_scale_per;
        end
        
    elseif c_type == 1
        disp('METTERE A POSTO')
%         % Gaussian gradient from center to periphery
%         side1 = size(atrial_tissue,1);
%         side2 = size(atrial_tissue,2);
%         [col, row] = meshgrid(1:side1, 1:side2);
%         sigma_x = x_r/2;
%         sigma_y = y_r*2/3*1/2;
%         gauss_bell = exp(-( ((col-x_c).^2/(2*sigma_x^2)) + (row-y_c).^2/(2*sigma_y^2))); %% gaussian function
%         gauss_bell = max(gauss_bell(:)) -gauss_bell ;
%         for i = 1: side1
%             for j = 1:side2
%                 if atrial_tissue(i,j) == 1                    
%                     g_scale_per = ones(size(Per_factor)) + (Per_factor-ones(size(Per_factor))) * gauss_bell(i,j);
%                     g_rand(i,j,:) = nomCond .* exp(sigma*randn(1,length(nomCond))) .* g_scale_per;
%                 else
%                     g_rand(i,j,:) = zeros(size(Per_factor));
%                 end
%             end
%         end
%         g_rand = reshape(g_rand, [], length(nomCond));
%         g_rand = g_rand(atrial_tissue == 1, :);
        
    elseif c_type == 2
        % sigmoidal distribution from center to peryphery
        side1 = size(atrial_tissue,1);
        side2 = size(atrial_tissue,2);
        [col, ~] = meshgrid(1:side1, 1:side2);
        x = (col - x_c);
        sig_boltz = 1./(1+exp(alpha* x + beta));
        sig_boltz(sep_high-sep_half_width:sep_high+sep_half_width-1, :) = sig_boltz(sep_high-sep_half_width:sep_high+sep_half_width-1, [4:end, end end end]);
        sig_boltz(sep_low-sep_half_width:sep_low+sep_half_width-1, :) = sig_boltz(sep_low-sep_half_width:sep_low+sep_half_width-1, [4:end, end end end]);

        
        for i = side1:-1:1
            for j = side2:-1:1
                if atrial_tissue(i,j) == 1                     
                    g_scale_per = ones(size(Per_factor)) + (Per_factor-ones(size(Per_factor))) * sig_boltz(i,j);
                    g_rand(i,j,:) = nomCond  .* g_scale_per; % .* exp(sigma*randn(1,length(nomCond)))
                else
                    g_rand(i,j,:) = zeros(size(Per_factor));
                end               
            end
        end 
        g_rand = reshape(g_rand, [], length(nomCond));
        g_rand = g_rand(atrial_tissue == 1, :);        
        
    end
    
end

end
