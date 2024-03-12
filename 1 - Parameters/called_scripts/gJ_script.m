
% San cells are connected to atrial cells with their value of gJ (1E9) ->
% quindi alla frontiera tra i due metto sempre 1e9, sia quando ho le
% cellule SAN, sia quando ho quelle atriali

%%% Check SAN/atrium frontier
gJ = atrial_tissue;
gJ_L = gJ - gJ(:,[1 1:end-1]);
gJ_U = gJ - gJ([1 1:end-1],:);
gJ_R = gJ - gJ(:,[2:end end]);
gJ_D = gJ - gJ([2:end end],:);

gJ_L = gJ_L(:);
gJ_U = gJ_U(:);
gJ_R = gJ_R(:);
gJ_D = gJ_D(:);

% Create 2D matrix [numcells, 4] con colonne L-U-R-D
gJ_matrix = zeros(side1*side2, 4);
gJ_matrix(gJ(:) == idx_san,:) = idx_san;
gJ_matrix(gJ(:) == idx_fat,:) = idx_fat;
gJ_matrix(gJ(:) == idx_fibro,:) = idx_fibro;

% manage border
gJ_matrix( abs(gJ_L) > 2, 1) = idx_fat;
gJ_matrix( abs(gJ_U) > 2, 2) = idx_fat;
gJ_matrix( abs(gJ_R) > 2, 3) = idx_fat;
gJ_matrix( abs(gJ_D) > 2, 4) = idx_fat;

%%% Choose uniform/gradient gJ distribution
if gJ_grad_flag == 0 % Uniform gJ distribution
    
    % Manage SAN/atrium interface -> SAME gJ values
    gJ_matrix( abs(gJ_R) == 1, 3) = idx_atr;
    gJ_matrix( abs(gJ_L) == 1, 1) = idx_atr;
    % Manage fibro/atrium interface -> SAME gJ values
    gJ_matrix( abs(gJ_R) == 3, 3) = idx_atr;
    gJ_matrix( abs(gJ_L) == 3, 1) = idx_atr;
    
    gJ_matrix(gJ_matrix == idx_atr) = G_atr;
    gJ_matrix(gJ_matrix == idx_san) = G_san;
    gJ_matrix(gJ_matrix == idx_fat) = 0;
    gJ_matrix(gJ_matrix == idx_fibro) = G_fib;
    
elseif gJ_grad_flag == 1 % Gradient
    if sep_type == 0
        %%% Gaussian Ggap distribution
        % sigma_x = 6;
        % sigma_y = y_r*2/3;
        % z = exp(-( ((col-x_c).^2/(2*sigma_x^2)) + (row-y_c).^2/(2*sigma_y^2))); %% gaussian function
        
        %%% Sigmoidal Ggap distribution
        % if length(grad_type) > 7 && strcmpi(grad_type(1:7), '-steep_')
        alfa_GJ = 1; % 1|10
        beta_GJ = -20; % -15|0
        % else
        %     alfa_sigm = 10; % 1|10
        %     beta_sigm = 0; % -15|0
        % end
        half_sep_x = x_c;
        half_sep_y = y_c-1;
        sigmoid_GJ =   G_atr + (G_san-G_atr)* 1./(1+exp((col - half_sep_x + beta_GJ ) /alfa_GJ ));
        sigmoid_GJ(sep_high-sep_half_width:sep_high+sep_half_width, :) = G_atr + (G_san-G_atr)* 1./(1+exp( (col(sep_high-sep_half_width:sep_high+sep_half_width, :) - half_sep_x + beta_GJ +3)/alfa_GJ ));
        sigmoid_GJ(sep_low-sep_half_width:sep_low+sep_half_width, :) = G_atr + (G_san-G_atr)* 1./(1+exp( (col(sep_low-sep_half_width:sep_low+sep_half_width, :) - half_sep_x + beta_GJ +3)/alfa_GJ ));
        
        % grad_shift = 88;
        % grad_steep = 5;
        % if strcmpi(grad_type, '-sweet_')
        %     half_sep_x = ;
        %     half_sep_y = y_c-1;
        %     beta_sigm = 0;
        %     sig_x = ((col - half_sep_x + beta_sigm) /(x_r*1.5)).^2 ;
        %     sig_y = ((row - half_sep_y) /(y_r*1.2)).^2 ;
        %     sigmoid_GJ = 1 ./ ( (1/G_atr) + (1.02e6./(1+exp(0.15*( sig_x + sig_y - 1 )))) );
        % else
        %     half_sep_x = grad_shift;
        %     half_sep_y = y_c-1;
        %     beta_sigm = 0;
        %     sig_x = ((col - half_sep_x + beta_sigm) /(x_r*1.5)).^2 ;
        %     sig_y = ((row - half_sep_y) /(y_r*1.2)).^2 ;
        %     sigmoid_GJ_ell = 1 ./ ( (1/G_atr) + ( (1/G_san) ./(1+exp(grad_steep*( sig_x + sig_y - 1 )))) );
        % end
    end
        
    % Manage SAN/atrium interface -> SAME gJ values
    gJ_matrix( abs(gJ_R) == 1, 3) = idx_san;
    gJ_matrix( abs(gJ_L) == 1, 1) = idx_san;
    % Manage fibro/atrium interface -> SAME gJ values
    gJ_matrix( abs(gJ_R) == 3, 3) = idx_atr;
    gJ_matrix( abs(gJ_L) == 3, 1) = idx_atr;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Manage fibro/SAN interface -> NO CURRENT EXCHANGE
    gJ_matrix( abs(gJ_R) == 2, 3) = idx_fat;
    gJ_matrix( abs(gJ_L) == 2, 1) = idx_fat;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % gJ_L is shifted by 1 so that gJ_D(i,j) == gJ_L(i+1, j)
    z_Lshifted = sigmoid_GJ(:, [1 1:end-1]);
    gJ_matrix(gJ_matrix(:, 1) == idx_san, 1) = gJ_matrix(gJ_matrix(:, 1) == idx_san, 1).* z_Lshifted(gJ_matrix(:, 1) == idx_san);
    for i = 2:4
        gJ_matrix(gJ_matrix(:, i) == idx_san, i) = gJ_matrix(gJ_matrix(:, i) == idx_san, i).*sigmoid_GJ(gJ_matrix(:, i) == idx_san);
    end
    
    gJ_matrix(gJ_matrix == idx_atr) = G_atr;
    gJ_matrix(gJ_matrix == idx_fat) = 0;
    gJ_matrix(gJ_matrix == idx_fibro) = G_fib;
end
