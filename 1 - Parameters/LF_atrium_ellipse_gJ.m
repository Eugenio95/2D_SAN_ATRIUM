close all
clear
clc

global sep_high sep_half_width sep_low

addpath called_scripts
addpath Seeds_hetero

adjust_size = 0; % input('SEP lunghe (0) o corte (10)?   '); % 0 default, 10 for short
SAN_cluster = 0;

nSEP = 5; % input('nSEP = ');
flag_species = 'human';
model_str    = 'Koivu'; %input('Choose model: "Koivu" or "mbs"  ', 's'); % 

sim_t           = 50; %% input('Sim duration in s = '); % seconds
sep_half_width  = 5; %% SEP opening width (total ? 2*sep_half_width+1)
gJ_grad_flag    = 3; %% input('Uniform (0) or SAN gradient (1) gap junction distribution or double sigmoidal (2) or SAN+atrium gradient (3)');
mosaic_flag     = 1; %% input('Mosaic model? (atrial cells in SEP: 0/1/2 NO/YES/square block in SEP)   ');
sep_type        = 0; %% input('Type of SEP: canal (0) or opening (1)? ');
sigma_SAN       = 0.2; %% input('Heterogeneity level in SAN: ');
cond_rand_flag  = 0; %% input('Random (0) or random+gradient (1) SAN maximal conductances distribution? ');
fib_flag        = 1; %% input('Insert SAN fibrosis (%) = ');
fib_connections = 1; %% 1: fibroblast have always Rgap = 1/G_san Ohm gap jucntions; 0: fibroblast follow gradient gap junctions
fat_flag        = 0; %% 1: add distributed fat (resistive barriers)
fat_flag_type   = 0; %% 0: All SAN; 1: only ellipse; 2: only SEPs

grad_type = [num2str(nSEP),'SEP', num2str(sep_half_width),'_']; % -steep_

%% General parameters
%%% Number of cells
side1 = 200;
side2 = 200;
%%% Sim duration
while isempty(sim_t)
    disp('Error: input duration of simulation')
    sim_t = input('Sim duration in s = '); % seconds
end

%%% Coupling resistance -> IN REALTà SONO CONDUTTANZE PERCHE' NEL CODICE CUDA ORA MOLTIPLICO PER GJ PER EVITARE DI DIVIDERE PER INFINITO
% REFERENCE VALUES:
% Amsaleg et al. usano 0.3114 (per cellula = 100um viene 32.2 KOhm) e 0.7060 S/m per atrio e CT intracellulare e
% 1.1186 e 2.5361 S/m extracellulare
G_atr = 1/1e6; %  %%% In Aslanidi et al 2009 usano 500nS (2 MOhm) in 3D rabbit
G_san = 1/1e7; % default = 1e9 Ohm -> too slow CV!!! Change to 1e8
G_fib = 1/1e8; % they express Cx45 as in SAN -> see lab notes

while isempty(sigma_SAN)
    disp('Error: input SAN heterogeneity level')
    sigma_SAN = input('Heterogeneity level in SAN: '); % seconds
end
sigma_atrio = 0;
sigma_fib = 0;

%%% Change Gmax randomization if hetero is present
if sigma_SAN > 0
    model_type = ['_Hetero', grad_type];
else
    model_type = ['_noHet', grad_type];
end

if strcmpi(flag_species, 'human')
    model_type = ['_HUMAN', model_type];
end

%%% Cell types
idx_atr = 0;
idx_san = 1;
idx_fibro = 3;
idx_fat = 9;

% Default seed
seed_number = 3; %input('Select seed number (1 to 5): ');
select_seed
model_type = [model_type, 'Seed', num2str(seed_number), '-'];

%% Geometry
geometry_script

%% Assign cell types
assign_cell_types_and_SEPs

%% Add resistive barriers (non-conducting fatty tissue)
if fat_flag == 1
    distributed_fat_amount = 0.1; % e.g. 0.1 = 10%
    add_distributed_fat
else
    fat_str = '';
end

%% gJ gradient
while isempty(gJ_grad_flag)
    disp('Error: input type gap junctional coupling')
    gJ_grad_flag = input('Uniform (0) or gradient (1) gap junction distribution? ');
end
if gJ_grad_flag == 0
    gJ_script
    grad_str = 'UNIF';
    disp('---> Uniform gap junctional distribution selected')
elseif gJ_grad_flag == 1
    gJ_script
    grad_str = 'GRAD';
    disp('---> Gradient gap junctional distribution selected')
elseif gJ_grad_flag == 2
    gJ_script_double_sigm
    grad_str = 'GRAD';
    disp('---> Separated SAN/atrial gradient gap junctional distribution selected')
elseif gJ_grad_flag == 3
    gJ_script_single_sigm
    grad_str = 'GRAD';
    disp('---> Separated SAN/atrial gradient gap junctional distribution selected')
end

%% Mosaic model -> insert atrial cells in SEPs
mosaic_model_script % always seed 1
select_seed % switch back to actual seed

%% Fibrosis (non-excitable but conductive fibroblasts)
%%% fibrosis_script
modified_fibrosis_script

if fib_flag ~= 0 && fib_connections == 1
    modify_fib_connections
end

%% Compute number of cells
n_san = length(find(atrial_tissue == idx_san));
n_atrial = length(find(atrial_tissue == idx_atr));
n_fat = length(find(atrial_tissue == idx_fat));

%% Display geometry parameters
fprintf(1, '\n')
disp(['numSAN = ', num2str(n_san)])
disp(['numATRIUM = ', num2str(n_atrial)])
disp(['numFAT = ', num2str(n_fat)])
disp(['numFIBRO = ', num2str(n_fibro)])

n_ratio = n_san/n_atrial * 100;
disp(['SAN/atrium ratio = ', num2str(round(n_ratio)), '%'])
fprintf(1, '\n')

%% Random conductances
block_ion = [1 1]; % (1) If e (2) ICaL 
randomize_conductance_script

%% Checks
%%% Geometry
aa = atrial_tissue;
aa(aa==1) = 25;
aa(aa==3) = 50;
aa(aa==9) = 75;
figure('Name','Geometria','NumberTitle','off');
imagesc(aa)
colormap([61 38 168; 17 190 185; 0 0 0;  255 255 19]/255) % 61 38 168

%%% Ggap distribution
figure
title_gJ = {'L', 'U', 'R' , 'D'}; 
for i = 1:4
    subplot(2, 2, i)
    imagesc(reshape(gJ_matrix(:,i)*1e9, side1, side2),'AlphaData',(reshape(gJ_matrix(:,i), side1, side2)) ~= 0)
    title(title_gJ{i})

    colormap(flipud(jet))
    cb = colorbar;
    title(cb, '[nS]')
    % caxis([1e-4 1e-1])
    % set(gca,'ColorScale','log')
end

%%% Ggap profile in SEPs
figure
aa = reshape(gJ_matrix(:,i), side1, side2);
plot(aa(100,:), 'LineWidth', 2)
set(gca, 'Yscale', 'log')
ylim([5e-10 2e-6])
title('Gap junctional conductances in SEPs')
xlabel('Cell # along x axis'), ylabel('G_{gj} [S]')
hold on
semilogy([112, 112], [1e-10 2e-6], 'k--', 'Linewidth', 1.5)

set(gca, 'FontWeight', 'bold', 'FontSize', 12)

%% Check % fibrosis in SEPs
SEP1 = atrial_tissue(35:45, 89:109);
SEP2 = atrial_tissue(65:75, 92:112);
SEP3 = atrial_tissue(95:105, 92:112);
SEP4 = atrial_tissue(125:135, 92:112);
SEP5 = atrial_tissue(155:165, 89:109);

fibro_in_SEP1 = ((numel(find(SEP1 == 3)))/(21*11))*100;
fibro_in_SEP2 = ((numel(find(SEP2 == 3)))/(21*11))*100;
fibro_in_SEP3 = ((numel(find(SEP3 == 3)))/(21*11))*100;
fibro_in_SEP4 = ((numel(find(SEP4 == 3)))/(21*11))*100;
fibro_in_SEP5 = ((numel(find(SEP5 == 3)))/(21*11))*100;

disp(['Average percent fibrosis in SEPs = ', num2str(round(mean([fibro_in_SEP1 fibro_in_SEP2 fibro_in_SEP3 fibro_in_SEP4 fibro_in_SEP5]), 1)), '%'])

%% Check fibrosis in SAN
cell_in_ellipse = atrial_tissue(find((col-x_c).^2./(x_r.^2) + (row-y_c).^2./(y_r.^2) <= 1));
SAN_in_ellipse = numel(find(cell_in_ellipse == 1));
fibro_in_ellipse = numel(find(cell_in_ellipse == 3));

disp(['Average percent fibrosis in SAN = ', num2str(round( fibro_in_ellipse/(SAN_in_ellipse+fibro_in_ellipse)*100,1 )), '%'])

%% Save sim parameters
% param_str = [model_type, num2str(sim_t), 's_', grad_str, mosaic_str, sep_str, '_F', num2str(fib_flag), fibro_str, fat_str, '_', cond_grad_string, datestr(date)];
param_str = [model_type, num2str(sim_t), 's_', grad_str, mosaic_str, 'F', num2str(fib_flag), fibro_str, fat_str, '_', cond_grad_string, datestr(date)];
sim_param = [datenum(date), gJ_grad_flag, sep_flag, fib_flag, side1, side2, sim_t, x_r, y_r, x_c, y_c, n_san, n_atrial, n_fat, n_fibro, round(n_ratio)]';

disp(param_str)

%%{
exp_fold = 'Sim_Param_Folder/'; % export folder
% Export gJ
writematrix(gJ_matrix, [exp_fold, 'gap_junc', param_str], 'Delimiter', 'space')
% Export geometry
writematrix(atrial_tissue(:), [exp_fold, 'atrial_geometry', param_str])
% Export rand cond
writematrix(rand_g(:), [exp_fold, 'rand_cond', param_str])
% Export parameters
writematrix(sim_param, [exp_fold, 'sim_param', param_str])
%%}

%% Set the true aspect ratios... COME SI FA? DIPENDE DALL'ORIENTAMENTO DELLE FIBRE
%{
figure(1)
daspect([1 3 1])

%% Plot fiber direction

X = 1:2:200;
Y = X;
U = zeros(size(atrial_tissue)/2);
V = ones(size(atrial_tissue)/2);

figure
imagesc(atrial_tissue)
hold on
quiver(X, Y, U, V, 'r')
%}

%% Check gJ stencil
% clc
% xx = 45;
% yy = 116;
% aa = sub2ind([200, 200], xx, yy);
% disp('        L        U         R         D')
% disp(gJ_matrix(aa, :))

%% Calcolo RAM occupata
%{
size_of_double = 8; % bytes
size_of_int = 4;

nCells = side1*side2;
nStates = 33;

size_Y0       = nCells*nStates*size_of_double;
size_Vidx     = nCells*size_of_int;
size_CellType = nCells*size_of_int;
size_Cm       = nCells*size_of_double;
size_Rand_g   = nCells*12*size_of_double;
size_gJ       = nCells*4*size_of_double;

total_RAM = size_Y0*2 + size_Vidx + size_CellType + size_Cm + size_Rand_g + size_gJ;
disp(['Occupied RAM at every timpestep = ', num2str(total_RAM/(2^20)), 'MB'])

%}
