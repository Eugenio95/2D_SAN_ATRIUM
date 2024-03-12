
% Atrial geometry
atrial_tissue = zeros(side1, side2);
[col, row] = meshgrid(1:side1, 1:side2);

% SAN ellipse geometry
x_c = floor(side1/2)-20; %default = -20;
y_c = floor(side2/2);
% x_r = floor(x_c/7);
y_r = y_c - floor(y_c/4); % default is /4
x_r = floor(y_r/6); % default is /6
ellipse = (col-x_c).^2./(x_r.^2) + (row-y_c).^2./(y_r.^2) <= 1;
% Fatty border geometry
if adjust_size == 0
    sep_stretch = 1.5;
else
    sep_stretch = 1; % default for human = 1.5, rabbit = 1
end
sep_right_h = x_c + 20*sep_stretch;    %lunghezza canali all'estremità
sep_right_m = x_c + floor(x_r+10)*sep_stretch;    %lunghezza canali medi
sep_right_c = x_c + floor(x_r+10)*sep_stretch;    %lunghezza canale centrale
sep_right_top = x_c + 16*sep_stretch;