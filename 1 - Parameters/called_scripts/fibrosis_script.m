
while isempty(fib_flag)
    disp('Error: input amount of fibrosis')
    fib_flag = input('Insert SAN fibrosis (%) = ');
end
if fib_flag > 0
    if sep_type == 0
        fib_SEP_flag = 3; %input('Fibroblasts inside SEPs? 0/1/2/3 (NO/YES/ONLY/GRADIENT)   ');
        while isempty(fib_SEP_flag)
            disp('Error: input if fibrosis required inside SEPs: 0/1/2 (NO/YES/ONLY)')
            fib_SEP_flag = input('Fibroblasts inside SEPs? 0/1/2 (NO/YES/ONLY)   ');
        end
    else
        fib_SEP_flag = 0;
    end
    
    if fib_SEP_flag == 0
        % only inside ellipse
        ind_ellipse = find((col-x_c).^2./(x_r.^2) + (row-y_c).^2./(y_r.^2) <= 1);
        ind_ellipse(ind_ellipse == 17701 | ind_ellipse == 22501) = [];
        fibro_str = 'A';
    elseif fib_SEP_flag == 1
        % also SEPs
        ind_ellipse = find(atrial_tissue == idx_san);
        fibro_str = 'B';
    elseif fib_SEP_flag == 2
        % Only in SEPs
        ind_ellipse = find(atrial_tissue == idx_san);
        ind_SAN = find((col-x_c).^2./(x_r.^2) + (row-y_c).^2./(y_r.^2) <= 1);
        ind_ellipse( ismember(ind_ellipse, ind_SAN) ) = [];
        
        fibro_str = 'C';
    elseif fib_SEP_flag == 3 % gradiente di fibrosi nelle SEP (Li+Fedorov 2023): Da 46% nel SAN a 28% nelle SEP (in media... quindi fino a 7?)     
        alfa_fib = 2; 
        beta_fib = -22; % in half-SEP, otherwise unbalanced % of fibrosis in SEP

        fib_in_SAN = 45.7; %45.7;
        fib_in_SEP = 9.7; %%%% 27.7*2-45.7, 1: average fibrosis in SEP
        
        half_sep_x = x_c;
        half_sep_y = y_c-1;
                
        col_fib = col';
        
        %%% Lineare
        if mosaic_flag == 0
            sigmoid_fib = [zeros(1, 68) -linspace(9.7, 45.7, 110-68)+(45.7+9.7), zeros(1, 90)];
        else
            sigmoid_fib = [zeros(1, 68) -linspace(9.7, 45.7, 105-68)+(45.7+9.7), zeros(1, 95)];
        end
        %%% Vera sigmoide
% % %         sigmoid_fib = fib_in_SEP + (fib_in_SAN-fib_in_SEP)* 1./(1+exp(( col_fib - half_sep_x + beta_fib) /alfa_fib ));        
% % %         sigmoid_fib(:, sep_high-sep_half_width:sep_high+sep_half_width) =  fib_in_SEP + (fib_in_SAN-fib_in_SEP)* 1./(1+exp(( col_fib(:, sep_low-sep_half_width:sep_low+sep_half_width) - half_sep_x + beta_fib +3)/alfa_fib ));
% % %         sigmoid_fib(:, sep_low-sep_half_width:sep_low+sep_half_width) =  fib_in_SEP + (fib_in_SAN-fib_in_SEP)* 1./(1+exp(( col_fib(:, sep_low-sep_half_width:sep_low+sep_half_width) - half_sep_x + beta_fib +3)/alfa_fib ));
        
        fib_tissue = atrial_tissue;
        for i = 1:side1
                fib_empty = randsample(side2, round(sigmoid_fib(:, i)*side2/100));
                fib_tissue(fib_empty, i) = 2;
        end
        ind_ellipse = find(fib_tissue == 2 & atrial_tissue == 1);
        
        fibro_str = 'D';
    end
    
    % Impose the same variability at every execution of the script
    ind_fib = ind_ellipse(randperm(length(ind_ellipse)));
    n_fibro = length(ind_fib);
    ind_fib = ind_fib(1:n_fibro);
    
    atrial_tissue( ind_fib ) =  idx_fibro;
    atrial_tissue = reshape(atrial_tissue, side1, side2);
else
    n_fibro = 0;
    fibro_str = '';
end
