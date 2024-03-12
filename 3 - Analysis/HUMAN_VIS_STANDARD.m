close all
clear
clc

addpath called_scripts/

addpath versione_HUMAN_20Sep2023/Results/
addpath versione_HUMAN_20Sep2023/Sim_param_folder/
% addpath versione_HUMAN-SAcRT_24Nov2023/Results/
% addpath versione_HUMAN-SAcRT_24Nov2023/Sim_param_folder/
addpath versione_HUMAN_20Sep2023/Results/ISO_baseline/


param_str = '50s_GRAD_M0_SEPc_F1D__06-Feb-2024';
type = '_HUMAN_Hetero5SEP5_Seed3-';
sacrt_str = '';
%%% OCCHIO: _SNRT va in sacrt_str
sim = '';
%%%
% name_str = [sacrt_str, type, sim, param_str];
name_str = [type, sim, param_str];

if strcmpi(sacrt_str, '_SNRT')
    movie_flag  = 1;
    SF_flag     = 0;
    all_SF_flag = 0;
    act_flag    = 1;
    sacrt_flag  = 1;
    save_flag   = 1;
    CV_flag     = 0;
    egm_flag    = 1;
    CL_flag     = 1;
else
    movie_flag  = 1;
    SF_flag     = 1;
    all_SF_flag = 1;
    act_flag    = 1;
    sacrt_flag  = 0;
    save_flag   = 1;
    CV_flag     = 1;
    egm_flag    = 1;
    CL_flag     = 1;
end

% load dei dati necessari per la visualizzazione
t = load(['t', name_str, '.txt']);
Vm = load([ 'Vm', name_str, '.txt']);
% I = load(['Igap', param_str, type,'.txt']);

param = load(['sim_param', [type, sim, param_str], '.txt']);
geom = load(['atrial_geometry', [type, sim, param_str], '.txt']);

disp(['Sim param date = ', datestr(param(1))])

AP_num_thresh = 1; % in last 2 s

nStates = 43;
idx_atr = 0;
idx_san = 1;
idx_fibro = 3;
idx_fat = 9;

if isempty(sim)
    result_folder_name = ['Results/', type(2:end-1), sim, sacrt_str];
else
    result_folder_name = ['Results/', type(2:end), sim(1:end-1), sacrt_str];
end
mkdir(result_folder_name)

%% Definizione di t (campionato) e Vm_mat
Vm_mat = reshape(Vm, length(Vm)/length(t), length(t));
clear Vm

idxAtrial = 1:param(5)*param(6);
idxSAN = find(geom == idx_san);

idxFAT = find(geom == idx_fat);
idxFIBRO = find(geom == idx_fibro);
idxAtrial([idxSAN; idxFAT; idxFIBRO]) = [];

Vm_SAN    = squeeze(Vm_mat(idxSAN, :));
Vm_atrium = squeeze(Vm_mat(idxAtrial, :));
% Vm_FAT    = squeeze(Vm_mat(idxFAT, :));
% Vm_fibro  = squeeze(Vm_mat(idxFIBRO, :));

% I_SAN    = squeeze(I_mat(idxSAN, :, :));
% I_FAT    = squeeze(I_mat(idxFAT, :, :));
% I_fibro  = squeeze(I_mat(idxFIBRO, :, :));
% I_atrium = squeeze(I_mat(idxAtrial, :, :));

%% Vm campioni presi su di una linea del nostro tessuto
geom_mat=reshape(geom,sqrt(length(geom)),sqrt(length(geom)));

[dim1,dim2]=size(geom);
dim=sqrt(dim1);
Vm_timeline=permute(reshape(Vm_mat,dim,dim,length(t)),[1,2,3]);

% evidenzio la riga presa in considerazione per i campioni
[r, c]=size(geom_mat);
l= 40;  % trovo la metà delle righe
linea=geom_mat(l,:);
figure('Name','Scanline','NumberTitle','off');
imagesc(geom_mat)
% axis on
hold on
plot([1 c],[l l],'g');
hold off

% si definiscono i valori campione per tipo di tessuto
% discriminazione delle cellule sulla linea l (_c sta per campione)

idxAtrial_c = 1:param(5);
geom_c=linea;
idx_SAN_c= find(geom_c == idx_san);
idx_FAT_c= find(geom_c == idx_fat);
idx_FIBRO_c=find(geom_c == idx_fibro);
idxAtrial_c([idx_SAN_c,idx_FAT_c,idx_FIBRO_c])=[];

Vm_c=Vm_timeline(l,:,:);
Vm_c=squeeze(Vm_c);

Vm_SAN_c=squeeze(Vm_c(idx_SAN_c, :));
Vm_FAT_c=squeeze(Vm_c(idx_FAT_c,:));
Vm_fibro_c=squeeze(Vm_c(idx_FIBRO_c,:));
Vm_atrium_c=squeeze(Vm_c(idxAtrial_c,:));

%% Visualizzazione delle cellule campione

% SAN
figure('Name','Vm nelle varie tipologie di cellule','NumberTitle','off')
hold on
subplot(221)
plot(t, Vm_SAN_c, 'g'), xlabel('Time [s]'),ylabel('Vm SAN [mV]'), title('SAN')
% ylim([-90 40])

% Atrio
subplot(222)
plot(t, Vm_atrium_c, 'b'), xlabel('Time [s]'),ylabel('Vm Atrio [mV]'), title('Atrium')
% ylim([-90 40])

% Fibro
if isempty(Vm_fibro_c) == 0
    subplot(223)
    plot(t, Vm_fibro_c,'r' ), xlabel('Time [s]'),ylabel('Vm fibro [mV]'), title('Fibroblasts')
    %     ylim([-90 40])
end

% FAT
subplot(224)
plot(t, Vm_FAT_c,'m' ), xlabel('Time [s]'),ylabel('Vm FAT [mV]'), title('Border')
% ylim([-90 40])

hold off
sample_name = [result_folder_name, '/sample', name_str,'.fig'];
saveas(gcf, sample_name)

clear Vm_c Vm_atrium_c Vm_SAN_c Vm_FAT_c Vm_fibro_c Vm_fibro Vm_FAT

%% Timeline della Vm
if movie_flag == 1
    if max(diff(t)) < 0.5e-3
        mov_undersamp = 100;
    else
        mov_undersamp = 10;
    end
    
    figure
    fig_folder_name = [result_folder_name, '/Movie'];
    mkdir(fig_folder_name)
    for i=1:mov_undersamp:length(t)
        
        imagesc(Vm_timeline(:,:,i))
        c=colorbar;
        c.Label.String='[mV]';
        caxis([-100 50])
        title(['t = ', num2str(t(i)), 's'])
        axis  image;
        drawnow
        
         fig_name = [result_folder_name, '/Movie','/mov', name_str, '_', num2str(i), '.jpg'];
         saveas(gcf, fig_name)
        
    end
end

%% EGM (both unipolar and bipolar)
if egm_flag == 1
    [egm_unipolar, egm_bipolar] = compute_egm(Vm_timeline, geom_mat);
    
    figure
    subplot(211)
    plot(t, egm_unipolar)
    xlabel('t [s]'), ylabel('Unipolar EGMs')
    legend('E1', 'E2')
    
    subplot(212)
    plot(t, egm_bipolar(1,:))
    xlabel('t [s]'), ylabel('Bipolar EGM 1-2')
    
    fig_egm= [result_folder_name, '/egm', name_str,'.fig'];
    saveas(gcf, fig_egm)
else
    egm_unipolar = [];
end

%% CL
clear Vm_timeline
if CL_flag == 1
    %%{
    noBeat_atrium = 0;
    noBeat_san = 0;
    
    % if strcmpi(sacrt_str, '_SAcRT')
    %     % Only last 2s
    %     t_last2s  = t(3000:end);
    %     Vm_atrium = Vm_atrium(:, 3000:end);
    %     Vm_SAN    = Vm_SAN(:, 3000:end);
    %     Vm_mat    = Vm_mat(:, 3000:end);
    % else
    % end
    
    disp('Extracting CL in atrium...')
    for i = size(Vm_atrium, 1):-1:1
        
        [dVdtmax, dVdtmax_pos{i}] = findpeaks(diff(Vm_atrium(i, :))./diff(t)', 'MinPeakHeight', 3e3, 'MinPeakDistance', 50 * 1e-3/min(diff(t)));
        
        [wrn, wrnID] = lastwarn;
        warning('off', wrnID)
        
        if length(dVdtmax_pos{i}) >= AP_num_thresh
            CL_atrium{i} = diff(t(dVdtmax_pos{i})) * 1e3;
        else
            CL_atrium{i} = nan;
            noBeat_atrium = noBeat_atrium + 1;
        end
    end
    
    %%%
    disp('Extracting CL in SAN...')
    for i = size(Vm_SAN, 1):-1:1
        [OS, OSpos{i}] = findpeaks(Vm_SAN(i, :), 'MinPeakHeight', 0, 'MinPeakDistance', 50 * 1e-3/min(diff(t)));
        
        [wrn, wrnID] = lastwarn;
        warning('off', wrnID)
        
        if length(OSpos{i}) >= AP_num_thresh
            CL_SAN{i} = diff(t(OSpos{i})) * 1e3;
            
        else
            CL_SAN{i} = nan;
            noBeat_san = noBeat_san +1;
        end
        
    end
    
    clc
    disp('CL extraction complete')
    
    disp(['Number of non beating SAN cells = ', num2str(noBeat_san), '/', num2str(size(Vm_SAN, 1))])
    disp(['Number of non beating atrial cells = ', num2str(noBeat_atrium), '/', num2str(size(Vm_atrium, 1))])
    
    if noBeat_san < 200 && noBeat_atrium < 1000
        
        CLsan_mean = round(mean( cell2mat(cellfun(@(x) x(end), CL_SAN, 'UniformOutput', false)), 'omitnan'));
        CLsan_std = round(std( cell2mat(cellfun(@(x) x(end), CL_SAN, 'UniformOutput', false)), 'omitnan'), 1);
        CLsan_min = round(min( cell2mat(cellfun(@(x) x(end), CL_SAN, 'UniformOutput', false)), [], 'omitnan'));
        disp(['CL in SAN = ', num2str(CLsan_mean), '±', num2str(CLsan_std), ' ms'])
        
        CLatrium_mean = round(mean( cell2mat(cellfun(@(x) x(end), CL_atrium, 'UniformOutput', false)), 'omitnan'));
        CLatrium_std = round(std( cell2mat(cellfun(@(x) x(end), CL_atrium, 'UniformOutput', false)) , 'omitnan'), 1);
        
        figure
        subplot(211)
        histogram(cell2mat(cellfun(@(x) x(end), CL_SAN, 'UniformOutput', false)))
        xlabel('CL [ms]')
        ylabel('Cell count')
        title('SAN')
        
        subplot(212)
        histogram(cell2mat(cellfun(@(x) x(end), CL_atrium, 'UniformOutput', false)))
        xlabel('CL [ms]')
        ylabel('Cell count')
        title('Atrium')
        
        hist_name = [result_folder_name, '/hist', name_str,'.fig'];
        saveas(gcf, hist_name)
        
        %%% SAN CL map
        CL_SAN_end_mat = cell2mat(cellfun(@(x) round(mean(x)), CL_SAN, 'UniformOutput', false));
        CL_atrium_end_mat = cell2mat(cellfun(@(x) round(mean(x)), CL_atrium, 'UniformOutput', false));
        CL_atrium_end_std = round(mean(cell2mat(cellfun(@(x) round(std(x), 1), CL_atrium, 'UniformOutput', false))),1);
        disp(['CL in atrium = ', num2str(mean(CLatrium_mean)), '±', num2str(CL_atrium_end_std), ' ms'])
        
        CL_map = geom_mat;
        CL_map(CL_map == 0) = CL_atrium_end_mat;
        CL_map(CL_map == 1) = CL_SAN_end_mat;
        % CL_map(CL_map ~= 1 | CL_map ~= 0) = nan;
        
        figure
        imagesc(CL_map)
        colormap([0 0 0; parula])
        c_CL = colorbar;
        c_CL.Label.String='[ms]';
        caxis([min(CL_SAN_end_mat)-1 max(CL_SAN_end_mat)])
        title('CL distribution')
        
        CLmap_name = [result_folder_name, '/CL', name_str,'.fig'];
        saveas(gcf, CLmap_name)
    end

else
    CL_SAN = [];
    CL_atrium = [];
    CLsan_mean = [];
    CLsan_std = [];
    CLatrium_mean = [];
    CLatrium_std = [];
end

clear Vm_atrium Vm_SAN

%% Safety factor
front_SEP_cells = find( geom_mat - geom_mat(:,[2:end end])  == 1) +200; % find the atrial cells of the frontier
tic
if SF_flag == 1 && noBeat_atrium < 1000
    
    %%% Carico file di gap junction
    gJ = reshape(load(['gap_junc', type, sim, param_str, '.txt']), size(geom_mat, 1), size(geom_mat, 2), 4);
    
    if all_SF_flag == 1
        front_SEP_cells = 1:size(Vm_mat, 1);
    else
        %%% Extract indexes of atrial cells at frontiers
        front_SEP_cells = find( geom_mat - geom_mat(:,[2:end end])  == 1) +200; % find the atrial cells of the frontier
        % aggiungo cellule a caso nell'atrio per testare SF
        %         front_SEP_cells = [front_SEP_cells', 201, 7940, 20000, 31840, 39840];
    end
    
    Igap_script
    
%     Igap = Igap *50; % serve solo nel SF Boyle & Vigmond
    
%         front2 = 1103; %21770, 22360; %1103; %8826; %22500
%         SFv2 = compute_SFv_tissue(front2, t, Vm_mat, Inet, Igap*50, geom_mat(:)) % versione Vigmond
%         SFz2 = compute_SFz_tissue(front2, t, Vm_mat, Inet-Igap(:, 1:end-1), Inet,  geom_mat) % versione Rudy
    %     SF2 = compute_SF_tissue(front2, t, Vm_mat, Inet, Igap,  geom_mat); % versione mista
    
    SF = compute_SFv_tissue(front_SEP_cells, t, Vm_mat, Inet, Igap*50, geom_mat(:)); % versione Vigmond
%     SFz = compute_SFz_tissue(front_SEP_cells, t, Vm_mat, Inet-Igap(:, 1:end-1), Inet,  geom_mat); % versione Rudy
    %     SF2 = compute_SF_tissue(front_SEP_cells, t, Vm_mat, Inet, Igap,  geom_mat); % versione mista
     
    %%% Display SF on tissue
    SF_map = geom_mat;
    SF_alpha = ones(size(SF_map));
    SF_alpha(SF_map == 1) = nan; % versione Vigmond solo atrio
    SF_alpha(SF_map == 9) = nan;
    SF_map(front_SEP_cells) = SF(:, 5);
    
    figure
    imagesc(SF_map, 'AlphaData', SF_alpha)
    colormap([0 0 0; cool])
    c_SF = colorbar;
    c_SF.Label.String='SF';
    caxis([0.99 3]) % caxis([min(SF)-1, max(SF)])
    %     caxis([min(SFv(:)), max(SFv(:))])
    set(gca,'color',[1 1 1])
    title('Safety factor')
    
    SFmap_name = [result_folder_name, '/SF', type, sacrt_str, param_str,'.fig'];
    saveas(gcf, SFmap_name)
    
    clear Igap Inet
    
else
    SF = [];
end
toc

%% CV nelle SEP
if CV_flag == 1
    for i = 0:10
        CV_sep1 = CV_in_SEP(t, Vm_mat, 35, 89, 67e-6, 20);
    end
    CV_sep2 = CV_in_SEP(t, Vm_mat, 70, 92, 67e-6, 20);
    for i = 0:6
        CV_sep3 = CV_in_SEP(t, Vm_mat, 100, 89, 67e-6, 20);
    end
    CV_sep4 = CV_in_SEP(t, Vm_mat, 130, 92, 67e-6, 20);
    CV_sep5 = CV_in_SEP(t, Vm_mat, 160, 89, 67e-6, 20);
    
    CV_san = CV_in_SEP_vert(t, Vm_mat, 100, 80, 67e-6, 30);
    CV_atr = CV_in_SEP(t, Vm_mat, 100, 135, 122.051e-6, 50);
else
    CV_sep1 = [];
    CV_sep2 = [];
   CV_sep3 = [];
   CV_sep4 = [];
   CV_sep5 = []; 
   CV_san = [];
   CV_atr = [];
end

%% SRRT
%%% prendo il primo AP in atrio dopo fine stimolazione (15.8 s)
%%% Così però nel caso di rientro è un po' un casino
%%% In realtà non è il SAcRT, ma potrei chiamarlo tipo
%%% "ATRIAL DRIVING RECOVERY TIME"
if sacrt_flag == 1
    figure, plot(t, Vm_mat(9900, :))
    
    t_SNRT_start = 41.5;
    t_last_stim = find(t > t_SNRT_start+0.2, 1, 'first');
%     t_SR = 16.3;
    
%     first_recovered_atrial_AP_idx = 99; 
    first_recovered_atrial_AP_idx = cell2mat(cellfun(@(x) find(x > t_last_stim, 1, 'first'), dVdtmax_pos(front_SEP_cells), 'UniformOutput', false));
    if all(first_recovered_atrial_AP_idx == first_recovered_atrial_AP_idx(1))
        first_recovered_atrial_AP = t(dVdtmax_pos{front_SEP_cells(1)}(first_recovered_atrial_AP_idx(1)));
%         first_recovered_atrial_AP = first_recovered_atrial_AP(find(first_recovered_atrial_AP > t_SR, 1, 'first'));
    else
        first_recovered_atrial_AP = t(dVdtmax_pos{1}(median(first_recovered_atrial_AP_idx)));
    end
    SNRT = (first_recovered_atrial_AP - t_SNRT_start)*1e3;
    disp(['SNRT = ', num2str(round(SNRT, 2)), ' s'])
else
    SNRT = [];
end

%% Activation map
%%{
front_SEP_cells = find( geom_mat - geom_mat(:,[2:end end])  == 1) +200; % find the atrial cells of the frontier

if act_flag == 1 && noBeat_atrium < 1000
    act_time_mat = get_activation_map(t, Vm_mat, geom_mat, front_SEP_cells);
    
    figure
    imagesc(act_time_mat)
    colormap([1 1 1; parula])
    c_AM = colorbar;
    c_AM.Label.String='[ms]';
    caxis([min(act_time_mat(:))- 10 max(act_time_mat(:))])
    set(gca,'color',[1 1 1])
    title('Activation map')
    
    activation_name = [result_folder_name, '/act_map', name_str,'.fig'];
    saveas(gcf, activation_name)   
end
%}

%% Save
if save_flag == 1
    save([result_folder_name, '/results_', name_str], 'CL_SAN', 'CL_atrium', 'CLsan_mean', 'CLsan_std', ...
        'CLatrium_mean', 'CLatrium_std', 'SF', 'SNRT', 'egm_unipolar',...
        'CV_sep1','CV_sep2','CV_sep3','CV_sep4','CV_sep5', 'CV_san','CV_atr')
    disp('Save complete')
end

%% LPM
%{
LPM_flag = 1;
if LPM_flag == 1
    [LPM_center, LPM_idx] = find_LPM_center(t, Vm_timeline, Vm_mat);
    disp(num2str(LPM_idx))
end
%}

%% CL in LPM identificati
%{
LPM1 = sub2ind([200, 200], 80, 106);
LPM2 = sub2ind([200, 200], 113, 106);
LPM3 = sub2ind([200, 200], 140, 106);

LPM1_san_idx = find(idxSAN == LPM1);
LPM2_san_idx = find(idxSAN == LPM2);
LPM3_san_idx = find(idxSAN == LPM3);

CL_LPM1 = mean(CL_SAN{LPM1_san_idx});
CL_LPM2 = mean(CL_SAN{LPM2_san_idx});
CL_LPM3 = mean(CL_SAN{LPM3_san_idx});
%}

%% Conduction velocity
%{
for i = size(Vm_mat, 1):-1:1

    [dVdtmax, dVdtmax_pos{i}] = findpeaks(diff(Vm_mat(i, :))./diff(t)', 'MinPeakHeight', 3e3, 'MinPeakDistance', 50);

    [wrn, wrnID] = lastwarn;
    warning('off', wrnID)
    
    if length(dVdtmax_pos{i}) >= AP_num_thresh
        CL_atrium{i} = diff(t(dVdtmax_pos{i})) * 1e3;
    else
        CL_atrium{i} = nan;
        noBeat_atrium = noBeat_atrium + 1;
    end
end

cell_start = sub2ind(size(geom_mat), 100, 135);
cell_end   = sub2ind(size(geom_mat), 100, 135+50);

if t(dVdtmax_pos{cell_end}(end)) > t(dVdtmax_pos{cell_start}(end))
    delta_t = t(dVdtmax_pos{cell_end}(end)) - t(dVdtmax_pos{cell_start}(end));
else
    delta_t = t(dVdtmax_pos{cell_end}(end)) - t(dVdtmax_pos{cell_start}(end-1));
end
delta_s = 50 * 100e-6; % 50 cell of 100 um length each

CV = round((delta_s / delta_t) * 100) ; % cm/s
disp(['CV = ', num2str(CV), ' cm/s'])
%}

%% Visualizzazione dati campione con Waterfall
%{
figure('Name', 'Vm per tipologia di cellula','NumberTitle','off')
hold on

% Waterfall cell SAN campione
subplot(221)
A=gradient(Vm_SAN_c);
waterfall(Vm_SAN_c,A);
colorbar, xlabel('Time [s]'),ylabel('n° of cells'), title('SAN'),zlim([-90 50])

% Waterfall cell FAT campione
subplot(222)
D=gradient(Vm_atrium_c);
waterfall(Vm_atrium_c,D);
colorbar, xlabel('Time [s]'),ylabel('n° of cells'), title('Atrium'),zlim([-90 50])

% Waterfall cell fibro campione
if isempty(Vm_fibro_c) == 0
    subplot(223)
    C=gradient(Vm_fibro_c);
    waterfall(Vm_fibro_c,C);
    colorbar, xlabel('Time [s]'), ylabel('n° of cells'),title('Fibroblasts'),zlim([-90 50])
end

% Waterfall cell atrium campione
subplot(224)
B=gradient(Vm_FAT_c);
waterfall(Vm_FAT_c,B);
colorbar, xlabel('Time [s]'),ylabel('n° of cells'), title('Border'),zlim([-90 50])
%}

%%
% prendo dei valori campione per poter semplificare la visualizzazione

%dimensioni dei vetori Vm_x
%dim SAN
% [row_s,col_s]=size(Vm_SAN);
% %dim FAT
% [row_fat,col_fat]=size(Vm_FAT);
% %dim fibro
% [row_fib,col_fib]=size(Vm_fibro);
% %dim atrium
% [row_a,col_a]=size(Vm_atrium);

% Vm_SAN_c=[Vm_SAN(1,:);Vm_SAN(floor(row_s/2),:);Vm_SAN(end,:)];
% Vm_FAT_c=[Vm_FAT(1,:);Vm_FAT(floor(row_fat/2),:);Vm_FAT(end,:)];
% Vm_fibro_c=[Vm_fibro(1,:);Vm_fibro(floor(row_fib/2),:);Vm_fibro(end,:)];
% Vm_atrium_c=[Vm_atrium(1,:);Vm_atrium(floor(row_a/2),:);Vm_atrium(end,:)];

%%
% close all
% figure
% plot(t, Vm_atrium, 'Color', [0.5,0.5,0.5]), xlabel('Time [s]'), ylabel('Vm [mV]')
% hold on
% plot(t, Vm_SAN, 'b')
% plot(t, Vm_fibro,'r')
% ylim([-90 40])
%
% figure('Name','Vm nelle varie tipologie di cellule','NumberTitle','off')
% hold on
% subplot(221)
% plot(t, Vm_SAN_c, 'g'), xlabel('Time [s]'),ylabel('Vm SAN [mV]'), title('SAN')
% ylim([-90 40])
%
% subplot(222)
% plot(t, Vm_atrium_c, 'b'), xlabel('Time [s]'),ylabel('Vm Atrio [mV]'), title('Atrium')
% ylim([-90 40])
%
% subplot(223)
% plot(t, Vm_FAT_c,'m' ), xlabel('Time [s]'),ylabel('Vm FAT [mV]'), title('Border')
% ylim([-90 40])
%
% subplot(224)
% plot(t, Vm_fibro_c,'r' ), xlabel('Time [s]'),ylabel('Vm fibro [mV]'), title('Fibroblasts')
% ylim([-90 40])

%% Plot Waterfall
% figure('Name', 'Vm per tipologia di cellula','NumberTitle','off')
% hold on
% %Waterfall 40 cell SAN
% subplot(221)
% A=gradient(Vm_SAN_c);
% waterfall(Vm_SAN_c,A);
% colorbar, xlabel('Time [s]'),ylabel('n° of cells'), title('SAN'),zlim([-90 50])
%
% %Waterfall 40 cell FAT
% subplot(222)
% B=gradient(Vm_FAT_c);
% waterfall(Vm_FAT_c,B);
% colorbar, xlabel('Time [s]'),ylabel('n° of cells'), title('Fibroblasts'),zlim([-90 50])
%
% %Waterfall 40 cell fibro
% subplot(223)
% C=gradient(Vm_fibro_c);
% waterfall(Vm_fibro_c,C);
% colorbar, xlabel('Time [s]'), ylabel('n° of cells'),title('Border'),zlim([-90 50])
%
%
% %Waterfall 40 cell atrium
% subplot(224)
% D=gradient(Vm_atrium_c);
% waterfall(Vm_atrium_c,D);
% colorbar, xlabel('Time [s]'),ylabel('n° of cells'), title('Atrium'),zlim([-90 50])
