
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

        fib_in_SAN = 45.7; %45.7; from Li et al. 2023 (Fedorov)
        fib_in_SEP = 9.7*1.5; %%%% 27.7*2-45.7 = 9.7, 1: average fibrosis in SEP
        fib_in_ATR = 0;
               
        fib_start_x = 89;
        
        %%% Lineare
        if adjust_size == 0
            sigmoid_fib = [zeros(1, fib_start_x-adjust_size) -linspace(fib_in_SEP, fib_in_SAN, 112-fib_start_x-adjust_size)+(fib_in_SAN+fib_in_SEP), -linspace(fib_in_ATR, fib_in_SEP, 5)+fib_in_SEP, zeros(1, 83+adjust_size)];
            sigmoid_fib(112:end-1) = sigmoid_fib(113:end);
        else
            sigmoid_fib = [zeros(1, fib_start_x-adjust_size) -linspace(fib_in_SEP, fib_in_SAN, 122-fib_start_x-adjust_size)+(fib_in_SAN+fib_in_SEP), -linspace(fib_in_ATR, fib_in_SEP, 5)+fib_in_SEP, zeros(1, 83+adjust_size)];
            sigmoid_fib(102:end-1) = sigmoid_fib(103:end);
        end
        
        sigmoid_fib = repmat(sigmoid_fib, 200, 1);        
        if SAN_cluster == 1              
            sigmoid_fib(:, 68:91) = 20;
        else                        
            sigmoid_fib(:, 68:91) = fib_in_SAN;
        end      
        sigmoid_fib([35:45 155:165], 1:end-3) = sigmoid_fib([35:45 155:165], 4:end);
        sigmoid_fib(:, 1:68) = 0;  
        
        %%% Aggiungo fibroblasti fuori dalle SEP
        atrial_tissue_extended = atrial_tissue;
        
        yy = [110 113, 113, 113, 110]-adjust_size;
        xx = [35, 65, 95, 125, 155];

        for i = 1:length(xx)
            atrial_tissue_extended(xx(i):xx(i)+2*sep_half_width, yy(i)-sep_half_width-1:yy(i)-1+5) = 1; % atrial cells in SEP have Rgap = 40 kOhm
        end
        %%%
        
        fib_tissue = atrial_tissue;
%         for i = 1:side1
%                 fib_empty = randsample(side2, round(sigmoid_fib(i, :)*side2/100));
%                 fib_tissue(fib_empty, i) = 2;
%         end
        
        fib_tissue(rand(200, 200)>(100-sigmoid_fib)/100 & atrial_tissue_extended == 1 & sigmoid_fib>0) = 2;
        ind_ellipse = find(fib_tissue == 2); % & atrial_tissue_extended == 1);
        
        fibro_str = 'D';
    end
    
    % Impose the same variability at every execution of the script
    ind_fib = ind_ellipse(randperm(length(ind_ellipse)));
    n_fibro = length(ind_fib);
    ind_fib = ind_fib(1:n_fibro);
    
    atrial_tissue( ind_fib ) =  idx_fibro;
    atrial_tissue = reshape(atrial_tissue, side1, side2);
        
    %%%%% CLUSTERS
    if SAN_cluster == 1
        im_clusters = imread('SAN_clusters.png');
        
        [length_pic, width_pic, depth_pic] = size(im_clusters);
        
        l1 = round(length_pic/200);
        l2 = round(width_pic/220);
        
        im_clusters_undersamp = im_clusters(2:l1:length_pic, 1:l2:width_pic, :);
        im_clusters_undersamp = im_clusters_undersamp(1:200, 1:200, :);
        im_clusters_undersamp(:, 1:end-5, :) = im_clusters_undersamp(:, 6:end,:);
        
        atrial_tissue(im_clusters_undersamp(:, :, 1) < 5) = idx_fibro;
        atrial_tissue(san_border == 1)  = idx_fat;
        
        n_fibro = length(find(atrial_tissue == 3));
    end
    
else
    n_fibro = 0;
    fibro_str = '';
end
