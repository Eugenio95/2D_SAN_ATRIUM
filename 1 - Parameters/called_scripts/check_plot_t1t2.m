
subplot(311)
plot(t, Vm_mat(cell_idx(i), :))

subplot(312)
plot(t(1:end-1), Inet(cell_idx(i), :))
hold on
plot(t(t2), Inet(cell_idx(i), t2), 'ko')
plot(t(t1), Inet(cell_idx(i), t1), 'r*')
title('Inet')
hold off

subplot(313)
plot(t, Igap(cell_idx(i), :))
hold on
plot(t(t2), Igap(cell_idx(i), t2), 'ko')
plot(t(t1), Igap(cell_idx(i), t1), 'r*')
title('Igap')
hold off