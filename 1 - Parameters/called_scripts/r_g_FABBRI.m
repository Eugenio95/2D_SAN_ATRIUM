function g_rand = r_g_FABBRI(sigma, num_cells, c_flag, c_type, alpha, beta, atrial_tissue, block, x_c)
    
global sep_high sep_half_width sep_low

% Geometry
Cm = 57;
L = 67;
R = 3.9;

% Nominal values
P_CaL = 0.4578;   % nanoA_per_millimolar (in i_CaL)
P_CaT = 0.04132;   % nanoA_per_millimolar (in i_CaT)
g_Kr = 0.00424;   % microS (in i_Kr)
K_NaCa = 3.343;   % nanoA (in i_NaCa)
i_NaK_max = 0.08105;   % nanoA (in i_NaK)
g_Na = 0.0223;   % microS (in i_Na)
g_Ks = 0.00065; % microS (in i_Ks)
g_f = 0.00427; % microS
g_to = 3.5e-3;   % microS (in i_to)
g_Kur = 0.1539e-3;   % microS (in i_Kur)
P_up_basal = 5.0;  % millimolar_per_second (in Ca_intracellular_fluxes)
g_KACh = 0.00345;   % microS (in i_KACh)

g_f = g_f * block(1);
P_CaL = P_CaL * block(2);

%% Nominal
% nomCond = [P_CaL, P_CaT, g_Kr, K_NaCa, i_NaK_max, g_Na, g_Ks, g_f, g_to, P_up_basal, g_KACh, Cm, L, R]; 
nomCond = [P_CaL, P_CaT, g_Kr, K_NaCa, i_NaK_max, g_Na, g_Ks, g_f, g_to, g_Kur, P_up_basal, g_KACh]; 

% Per_factor = [8, 3.2, 5, 10, 6.5, 10, 5, 1, 2.03, 2.03, 2.03]/2.03; %% CUSTOM -> gNa increase from (Li, et al. 2020)
Per_factor = [11.4, 3.2, 5, 10, 6.5, 20, 5, 12.6, 2.03, 2.03, 2.03]/2.03; %% From Ly & Weinberg (2018)
% Per_factor = [11.4, 3.2, 5, 10, 6.5, 20, 5, 12.6, 1, 1, 1 2.03, 1.23, 1.5]; %% From Ly & Weinberg (2018)
% where = 1 no data about that, but ony random hetero

%% Randomization
if c_flag == 0
    
    % Random hetero, no influence on space
    for i = num_cells:-1:1
        g_rand(i,:) = nomCond .* exp( sigma * randn(1,length(nomCond)) );
    end
%     g_rand = [g_rand, zeros(size(g_rand, 1), 1)];
    
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
        sig_boltz = 1./(1+exp((col - x_c + beta)/alpha));
        sig_boltz(sep_high-sep_half_width:sep_high+sep_half_width, :) = sig_boltz(sep_high-sep_half_width:sep_high+sep_half_width, [4:end, end end end]);
        sig_boltz(sep_low-sep_half_width:sep_low+sep_half_width, :) = sig_boltz(sep_low-sep_half_width:sep_low+sep_half_width, [4:end, end end end]);

        
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
