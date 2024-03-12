while isempty(cond_rand_flag)
    disp('Error: input of heterogeneity')
    cond_rand_flag = input('Random (0) or random+gradient (1) SAN maximal conductances distribution? ');
end
if cond_rand_flag == 0
    disp('---> Random Gmax distribution selected')
    cond_grad_type = [];
    %     if sigma_SAN > 0
    %         cond_grad_string = ''; % 'RAND_';
    %     else
    cond_grad_string = '';
    %     end
else
    % Heterogeneity (random Gmax) with peripheral phenotype (gradient)
    disp('---> Random+gradient Gmax distribution selected')
    disp('')
    disp('Linear superior/inferior gradient (0)')
    disp('or Exponential center/periphery gradient(1)')
    cond_grad_type = 2; %input('or Sigmoidal center/peiphery (2)? ');
    switch cond_grad_type
        case 0
            cond_grad_string = 'RAND+SupInf_';
        case 1
            cond_grad_string = 'RAND+EXP_';
        case 2
            cond_grad_string = 'RAND+SIGM_';
    end
end

alfa_sigm_pSAN = -1;
beta_sigm_pSAN = -22;

% block_ion = 1; %input('block = ');

if strcmpi(flag_species, 'rabbit')
    g_SAN = r_g_SEV(sigma_SAN, n_san, cond_rand_flag, cond_grad_type, alfa_sigm_pSAN, beta_sigm_pSAN, atrial_tissue, block_ion, x_c, sep_high, sep_half_width, sep_low);
    g_ATRIUM = r_g_LIND(sigma_atrio, n_atrial, g_SAN);
elseif strcmpi(flag_species, 'human')
    g_SAN = r_g_FABBRI(sigma_SAN, n_san, cond_rand_flag, cond_grad_type, alfa_sigm_pSAN, beta_sigm_pSAN, atrial_tissue, block_ion, x_c);
    if strcmpi(model_str, 'koivu')
        g_ATRIUM = r_g_KOIVU(sigma_atrio, n_atrial, g_SAN);
    elseif strcmpi(model_str, 'mbs')
        g_ATRIUM = r_g_MBS(sigma_atrio, n_atrial, g_SAN);
    end
end

rand_g = zeros(side1*side2, max(size(g_SAN, 2), size(g_ATRIUM, 2)));
rand_g(atrial_tissue(:) == idx_san,:) = [g_SAN, zeros(size(g_SAN, 1), size(rand_g, 2)-size(g_SAN, 2))];
rand_g(atrial_tissue(:) == idx_atr,:) = [g_ATRIUM, zeros(size(g_SAN, 1), size(rand_g, 2)-size(g_ATRIUM, 2))];

if n_fibro > 0
    g_MOR = r_g_MOR(sigma_fib, n_fibro);
    rand_g(atrial_tissue(:) == idx_fibro,:) = g_MOR;
end

rand_g = rand_g';
