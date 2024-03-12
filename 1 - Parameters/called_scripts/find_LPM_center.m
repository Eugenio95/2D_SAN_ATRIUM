function [LPM_center, LPM_sub] = find_LPM_center(t, Vm_timeline, Vm_mat)

figure
imagesc(Vm_timeline(:,:,1000))
c=colorbar;
c.Label.String='[mV]';
caxis([-100 50])
title(['t =', num2str(t(1000)), 's'])
axis  image;

[LPM_center(2), LPM_center(1)] = ginput(1);
LPM_center = round(LPM_center);
square_side = 3;
[A,B] = meshgrid(LPM_center(1)-square_side:LPM_center(1)+square_side, LPM_center(2)-square_side:LPM_center(2)+square_side);
c = cat(2,A',B');
d = reshape(c,[],2);
LPM_square = sub2ind([200, 200], d(:, 1), d(:, 2));

LPM_traces = reshape(Vm_mat(LPM_square, :), [], 2000);
plot(t, LPM_traces)

% Find first peak
for i = 1:size(LPM_traces, 1)
    [~, OSpos_LPM] = findpeaks(LPM_traces(i, :), 'MinPeakHeight', 0, 'MinPeakDistance', 100);
    if ~isempty(OSpos_LPM)
        firstOSpos(i) = OSpos_LPM(1);
    end
end

% Find peaks after first standard deviation low point (more similar curves)
[~, std_pos] = findpeaks(-std(LPM_traces), 'MinPeakProminence', 5, 'MinpeakDistance', 200);

firstOSpos(firstOSpos < std_pos(1)) = nan;
LPM_idx = LPM_square(find(firstOSpos == min(firstOSpos), 1));
[LPM_sub(1), LPM_sub(2)] = ind2sub([200, 200], LPM_idx);

hold on
plot(t, LPM_traces(find(firstOSpos == min(firstOSpos), 1), :), 'k-*')

end
