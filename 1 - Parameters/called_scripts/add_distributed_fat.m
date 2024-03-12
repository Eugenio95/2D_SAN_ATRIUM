
if fat_flag_type == 0 % Uniformly distributed in all SAN (ellipse + SEPs)
    ind_ellipse = find(atrial_tissue == idx_san);
    
    distributed_fat_prob = rand(1, length(ind_ellipse));
    distributed_fat_vector = ind_ellipse(distributed_fat_prob > (1-distributed_fat_amount));
    atrial_tissue(distributed_fat_vector) = 9;
    fat_str = 'FAT-A';
    
elseif fat_flag_type == 1 % Uniformly distributed only in SAN ellipse    
    ind_ellipse = find((col-x_c).^2./(x_r.^2) + (row-y_c).^2./(y_r.^2) <= 1);
    ind_ellipse(ind_ellipse == 17701 | ind_ellipse == 22501) = [];
    
    distributed_fat_prob = rand(1, length(ind_ellipse));
    distributed_fat_vector = ind_ellipse(distributed_fat_prob > (1-distributed_fat_amount));
    atrial_tissue(distributed_fat_vector) = 9;
    fat_str = 'FAT-B';
    
elseif fat_flag_type == 2 % Uniformly distributed only in SEPs    
    ind_ellipse = find(atrial_tissue == idx_san);
    ind_SAN = find((col-x_c).^2./(x_r.^2) + (row-y_c).^2./(y_r.^2) <= 1);
    ind_ellipse( ismember(ind_ellipse, ind_SAN) ) = [];
    
    
    distributed_fat_prob = rand(1, length(ind_ellipse));
    distributed_fat_vector = ind_ellipse(distributed_fat_prob > (1-distributed_fat_amount));
    atrial_tissue(distributed_fat_vector) = 9;       
    fat_str = 'FAT-C';

end

