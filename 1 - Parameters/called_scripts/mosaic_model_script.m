
while isempty(mosaic_flag)
    disp('Error: input if atrial cells are required in SEP: 0/1 NO/YES)')
    mosaic_flag = input('Mosaic model? (atrial cells in SEP: 0/1 NO/YES)   ');
end

%}     

load('Seed1.mat');
rng(s)

switch  mosaic_flag
    case 0
        mosaic_str = '_M0_';
        
    case 1
        alfa_MO = 2; %%% 2 o 1e-10 %%% 3 per sim completa
        beta_MO = -29;
        sigmoid_MO = 1./(1+exp( (col - x_c  + beta_MO)/ alfa_MO));
        sigmoid_MO(sep_high-sep_half_width:sep_high+sep_half_width, :) = 1./(1+exp( (col(sep_high-sep_half_width:sep_high+sep_half_width, :) - x_c + beta_MO +3)/alfa_MO ));
        sigmoid_MO(sep_low-sep_half_width:sep_low+sep_half_width, :) = 1./(1+exp( (col(sep_low-sep_half_width:sep_low+sep_half_width, :) - x_c + beta_MO +3)/alfa_MO ));
        
        if alfa_MO < 1e-3
            sigmoid_MO = floor(sigmoid_MO);
        end
        
        atrial_tissue(atrial_tissue == idx_san) = round( sigmoid_MO(atrial_tissue == idx_san) + 1* rand(size(atrial_tissue(atrial_tissue == idx_san))) -0.5 );
        
        mosaic_str = ['_M1',num2str(seed_number),'_'];
        
    case 2
        alfa_MO = 1e-10; %1e-10; %%% 2 o 1e-10
        beta_MO = -29;
        sigmoid_MO = 1./(1+exp( (col - x_c  + beta_MO)/ alfa_MO));
        sigmoid_MO(sep_high-sep_half_width:sep_high+sep_half_width, :) = 1./(1+exp( (col(sep_high-sep_half_width:sep_high+sep_half_width, :) - x_c + beta_MO +3)/alfa_MO ));
        sigmoid_MO(sep_low-sep_half_width:sep_low+sep_half_width, :) = 1./(1+exp( (col(sep_low-sep_half_width:sep_low+sep_half_width, :) - x_c + beta_MO +3)/alfa_MO ));
        
        if alfa_MO < 1e-3
            sigmoid_MO = floor(sigmoid_MO);
        end
        
        atrial_tissue(atrial_tissue == idx_san) = round( sigmoid_MO(atrial_tissue == idx_san) + 1* rand(size(atrial_tissue(atrial_tissue == idx_san))) -0.5 );
        
        mosaic_str = '_M2_';
        
end
