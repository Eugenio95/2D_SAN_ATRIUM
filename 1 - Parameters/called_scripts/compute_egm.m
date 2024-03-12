function [egm_unipolar, egm_bipolar] = compute_egm(Vm, geom_mat)

%%% Input parameters

%%% 2 coordinates for bipolar acquisition: the electrodes are 30 cells =
%%% 2 mm away from each other (Li et al. 2017 "Redundant ...")
xc = [60, 60, 60, 60]; % x coordinate of catheter tip
yc = [50, 70, 120, 140]-20; % y coordinate of catheter tip

[node(:, 1), node(:, 2)] = ind2sub([size(Vm, 1), size(Vm, 2)], 1:size(Vm,1)*size(Vm,2)); % coordinate di ogni cellula?
distance_surface = .25; % distance of catheter from atrial tissue in mm
tmax = size(Vm, 3); % # of samples in time

% Queste sono le distanze tra gli elettrodi? Se ne ho uno solo?
% Oppure sono i passi di discretizzazione della mesh?
setting.dx = 3; % mm?
setting.dy = 3;
%%% QUESTI SONO I COEFFICIENTI DI DIFFUSIONE in cm^2/s [0.001 in Shillieto et al 2016]
setting.Dx = 3;
setting.Dy = 3;

%%% Plot electrodes on atrial tissue geometry
figure
imagesc(geom_mat)
hold on
plot(xc,yc,'or','MarkerSize',8,'MarkerFaceColor','r');

%%% Computation
n_ele = length(xc); %number of electrodes
egm_unipolar = zeros(n_ele,tmax);
%registrazione monopolare si considera un elettrodo attivo, posto su una regione
%in cui ci sono variazione nel potenziale d’interesse, ed un altro posto su una regione che dovrebbe essere a potenziale 0.

for i = 1:n_ele
    grid_distance = sqrt((xc(i)-node(:,1)).^2+(yc(i)-node(:,2)).^2+distance_surface.^2);
    %   grid_distance = sqrt((xc(i)-node(:,1)).^2+(yc(i)-node(:,2)).^2)+eps;
    for t = 1:tmax
        % Conf Proc IEEE Eng Med Biol Soc. 2016 August ; 2016: 2741–2744. doi:10.1109/EMBC.2016.7591297.
        % Catheter Simulator Software Tool to Generate Electrograms of Any Multi-polar Diagnostic Catheter from 3D Atrial Tissue
        % Kristina E. Shillieto*, Prasanth Ganesan, M.S.*, Anthony J. Salmin*, Elizabeth M. Cherry, Ph.D.*, Arkady M. Pertsov, Ph.D.+, and Behnaz Ghoraani, Ph.D.*,†
        
        X = zeros(size(Vm,1)+2);
        X(2:end-1,2:end-1) = Vm(:,:,t);
        Vxx = diff(X,2,2)./(setting.dx^2);
        Vyy = diff(X,2,1)./(setting.dy^2);
        grad = 0.5 * setting.Dx * Vxx(2:end-1,:)  + 0.5 * setting.Dy * Vyy(:,2:end-1); %flip(
        egm_unipolar(i,t) = sum(grad(:)./(grid_distance.^2));
                       
        % Add noise? -> SNR = 50 dB
        
        if mod(t, 1000) == 0
            clc
            disp(['Computing EGM for electrode ', num2str(i), '/', num2str(n_ele), ': ', num2str(round((t)/tmax*100)), '%'])
        end
    end
end

% registrazione bipolare
egm_bipolar = zeros(n_ele,tmax);
for i = 1:n_ele-1
    egm_bipolar(i,:) = egm_unipolar(i,:) - egm_unipolar(i+1,:);
end

clc
disp('EGM computation complete')
    
end

