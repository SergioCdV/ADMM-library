
fig = openfig('..\..\..\Escritorio\Research\Publications\Press\RdVADMM_2023_EUCASS\Images\CW\L2\L2_D_cost.fig', 'reuse');
lines = findobj(fig, 'Type', 'Line');

% Si quieres ver cuántos hay
disp(length(lines));

% Extrae los datos de la primera curva, por ejemplo
x = get(lines(1), 'XData');
y = get(lines(1), 'YData');

% Dibuja en escala log-log
figure; 
loglog(x, y);
grid on;
xlabel('Iteration $i$');
ylabel('$\Delta V_T$ [m/s]');
