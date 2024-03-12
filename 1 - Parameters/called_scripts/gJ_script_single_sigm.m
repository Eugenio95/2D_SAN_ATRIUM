
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

% Manage border
gJ_matrix( abs(gJ_L) > 2, 1) = idx_fat;
gJ_matrix( abs(gJ_U) > 2, 2) = idx_fat;
gJ_matrix( abs(gJ_R) > 2, 3) = idx_fat;
gJ_matrix( abs(gJ_D) > 2, 4) = idx_fat;

n = 8*sep_half_width+1;
shift_y = -sep_half_width+8;

if sep_type == 0
    %%% Sigmoidal Ggap distribution
    % if length(grad_type) > 7 && strcmpi(grad_type(1:7), '-steep_')
    alfa_GJ = ((n+shift_y)/(41))*2; % 1
    beta_GJ = -40 + shift_y; % -25; default = -40 for human
    
    half_sep_x = x_c;
    half_sep_y = y_c-1;
    sigmoid_GJ =   G_atr + (G_san-G_atr)* 1./(1+exp((col - half_sep_x + beta_GJ ) /alfa_GJ ));
    sigmoid_GJ(sep_high-sep_half_width:sep_high+sep_half_width, :) = G_atr + (G_san-G_atr)* 1./(1+exp( (col(sep_high-sep_half_width:sep_high+sep_half_width, :) - half_sep_x + beta_GJ +3 +adjust_size)/alfa_GJ ));
    sigmoid_GJ(sep_low-sep_half_width:sep_low+sep_half_width, :) = G_atr + (G_san-G_atr)* 1./(1+exp( (col(sep_low-sep_half_width:sep_low+sep_half_width, :) - half_sep_x + beta_GJ +3+adjust_size)/alfa_GJ ));
end

% Manage SAN/atrium interface -> SAME gJ values
gJ_matrix( abs(gJ_R) == 1, 3) = idx_san;
gJ_matrix( abs(gJ_L) == 1, 1) = idx_san;
% Manage fibro/atrium interface -> SAME gJ values
gJ_matrix( abs(gJ_R) == 3, 3) = idx_atr;
gJ_matrix( abs(gJ_L) == 3, 1) = idx_atr;
% Manage fat/SAN interface -> NO CURRENT EXCHANGE
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

yy = [110 113, 113, 113, 110] - adjust_size;
xx = [35, 65, 95, 125, 155];

beta_x = ceil(n/2);
r_x = floor(beta_x/2);

y0 = sigmoid_GJ(1, 112:132)*1e6; % for MegaOhm

[rx,ry]= meshgrid( -n/2:n/2-1 , -n/2:n/2);
r = sqrt((0.5*rx).^2+(ry).^2);
vq = interp1( 0:(length(y0)-1) , y0 , r(:) , 'linear' , 0 );
W=r; W(:)=vq;
W = fliplr(W(1:round(end/2), :)');
W(W == 0) = 1;
W([1 round(end/2)], :) = [];
sigmoid_RA = repmat(W, 1, 1, 4); %/(max(sigmoid_RA(:)));

%%% Adjust Rgap in height out of the SEPs
%%% Atrial cells above half SEP height connect to cell below with Rgap of 
%%% the cell below and vicecersa
sigmoid_RA(1:round(end/2-1),  :, 4) = sigmoid_RA(2:(round(end/2)),  :, 4); % upper half down
sigmoid_RA(round(end/2+1):end,  :, 2) = sigmoid_RA(round(end/2):end-1,  :, 2); % lower half up

sigmoid_RA(end, :, 4) = 1;  % lower border of FUNGO as atrial cells
sigmoid_RA(1, :, 2) = 1;  % lower border of FUNGO as atrial cells

for j = 1:size(sigmoid_RA, 3)
    for i = 1:length(xx)
        % manage Up and Down gJ_matrices
        if j == 1
            shift_left = -1;
            start = 0;
            stop = 0;
        elseif j == 2
            shift_left = 0;
            start = 1;
            stop = 0;
        elseif j == 3
            shift_left = 0;
            start = 0;
            stop = 0;
        elseif j == 4
            shift_left = 0;
            start = 0;
            stop = -1;
        end
        
        idx_new_x = (xx(i)-floor((beta_x + r_x)/2-1):xx(i)-floor((beta_x + r_x)/2)+(2*beta_x-3)); %  + shift_y;
        idx_new_y = yy(i)-shift_left:yy(i)-shift_left+(beta_x-1);

        idx_new   = sub2ind([200, 200], repmat(idx_new_x, idx_new_y(end)-idx_new_y(1)+1, 1), repmat(idx_new_y', 1, length(idx_new_x)));
        sig = sigmoid_RA(:, :, j)'; %reshape(permute(sigmoid_RA, [2, 1, 3]), size(sigmoid_RA,1)*size(sigmoid_RA,2), 4);
        
        if sep_half_width > 6
            switch sep_half_width
                case 7
                    ovlap = 6;
                case 8
                    ovlap = 11;
                case 9
                    ovlap = 16;
                case 10
                    ovlap = 20;
            end
            if i == 1
                sig( 1:3, 1:ovlap) = 1;
            elseif i == 5
                sig( 1:3, (end-ovlap):end) = 1;
            end
            
            sig = fliplr(sig);
        end
        
        gJ_matrix( idx_new(:), j) = ...
            gJ_matrix(  idx_new(:), j) .* sig(:); % atrial cells outside of SEPs have circular sigmoidal gradient from 40 kOhm to 1 kOhm        
    
    end
end
