
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

if sep_type == 0
    %%% Sigmoidal Ggap distribution
    % if length(grad_type) > 7 && strcmpi(grad_type(1:7), '-steep_')
    alfa_GJ = 1; % 1
    beta_GJ = -25; % -25
    
    half_sep_x = x_c;
    sigmoid_GJ =   G_CT + (G_san-G_CT)* 1./(1+exp((col - half_sep_x + beta_GJ ) /alfa_GJ ));
    sigmoid_GJ(sep_high-sep_half_width:sep_high+sep_half_width, :) = G_CT + (G_san-G_CT)* 1./(1+exp( (col(sep_high-sep_half_width:sep_high+sep_half_width, :) - half_sep_x + beta_GJ +3)/alfa_GJ ));
    sigmoid_GJ(sep_low-sep_half_width:sep_low+sep_half_width, :) = G_CT + (G_san-G_CT)* 1./(1+exp( (col(sep_low-sep_half_width:sep_low+sep_half_width, :) - half_sep_x + beta_GJ +3)/alfa_GJ ));
    
    % Manage SAN/atrium interface -> SAME gJ values
    gJ_matrix( abs(gJ_R) == 1, 3) = idx_san;
    gJ_matrix( abs(gJ_L) == 1, 1) = idx_san;
    % Manage fibro/atrium interface -> SAME gJ values
    gJ_matrix( abs(gJ_R) == 3, 3) = idx_atr;
    gJ_matrix( abs(gJ_L) == 3, 1) = idx_atr;
    
    % Manage fibro/SAN interface -> NO CURRENT EXCHANGE
    gJ_matrix( abs(gJ_R) == 2, 3) = idx_fat;
    gJ_matrix( abs(gJ_L) == 2, 1) = idx_fat;
    
    % gJ_L is shifted by 1 so that gJ_D(i,j) == gJ_L(i+1, j)
    z_Lshifted = sigmoid_GJ(:, [1 1:end-1]);
    gJ_matrix(gJ_matrix(:, 1) == idx_san, 1) = gJ_matrix(gJ_matrix(:, 1) == idx_san, 1).* z_Lshifted(gJ_matrix(:, 1) == idx_san);
    for i = 2:4
        gJ_matrix(gJ_matrix(:, i) == idx_san, i) = gJ_matrix(gJ_matrix(:, i) == idx_san, i).*sigmoid_GJ(gJ_matrix(:, i) == idx_san);
    end
    
    gJ_matrix(gJ_matrix == idx_atr) = G_atr;
    gJ_matrix(gJ_matrix == idx_fat) = 0;
    gJ_matrix(gJ_matrix == idx_fibro) = G_fib;
    
    %%% FINO A QUI UGUALE A gJ_script
    %%% In più (FUNGO ATTORNO A SEPs):
    yy = [111 114, 114, 114, 111];
    xx = [35, 65, 95, 125, 155];
    
    alfa = -alfa_GJ;
    beta_x = 21;
    beta_y = 0;
    
    r_x = floor(beta_x/2);
    r_y = r_x;
    
    [col_RA, row_RA] = meshgrid(1:beta_x, 1:(2*beta_x-1));
    sigmoid_RA = 1- ((1-G_CT*1e3)*(( 1./ 1+ exp( ( ((col_RA - beta_y)/ r_y).^2 + ((row_RA - beta_x)/ (r_y)).^2 ) /alfa) )-1));
    sigmoid_RA = repmat(sigmoid_RA, 1, 1, 4);
    
    sigmoid_RA(:, :, 2) = sigmoid_RA([1 1:end-1],  :, 2);
    sigmoid_RA(1, :, 2) = 1;    % upper border of FUNGO as atrial cells
    sigmoid_RA(end, :, 4) = 1;  % lower border of FUNGO as atrial cells
    
    gJ_matrix = reshape(gJ_matrix, side1, side2, 4);
    for j = 1:size(gJ_matrix, 3)
        for i = 1:length(xx)
            % manage Up and Down gJ_matrices
            if j == 1
                shift_left_SEP = 0;
                shift_left = 0;
                start = 0;
                stop = 0;
            elseif j == 2
                shift_left_SEP = 1;
                shift_left = 1;
                start = 1;
                stop = 0;
            elseif j == 3
                shift_left_SEP = 1;
                shift_left = 1;
                start = 0;
                stop = 0;
            elseif j == 4
                shift_left_SEP = 1;
                shift_left = 1;
                start = 0;
                stop = -1;
            end
            
            if mosaic_flag == 2
                gJ_matrix(xx(i)+start:xx(i)+stop+2*sep_half_width, yy(i)-sep_half_width-shift_left_SEP:yy(i)-1-shift_left, j) = G_CT; % atrial cells in SEP have Rgap = 40 kOhm
            end
            gJ_matrix(xx(i)-floor((beta_x + r_x)/2):xx(i)-floor((beta_x + r_x)/2)+(2*beta_x-2), yy(i)-shift_left:yy(i)-shift_left+(beta_x-1), j) = ...
                gJ_matrix(xx(i)-floor((beta_x + r_x)/2):xx(i)-floor((beta_x + r_x)/2)+(2*beta_x-2), yy(i)-shift_left:yy(i)-shift_left+(beta_x-1), j) .* sigmoid_RA(:, :, j); % atrial cells outside of SEPs have circular sigmoidal gradient from 40 kOhm to 1 kOhm
            
        end
    end    
    
end

gJ_matrix = reshape(gJ_matrix, side1 * side2, 4);
