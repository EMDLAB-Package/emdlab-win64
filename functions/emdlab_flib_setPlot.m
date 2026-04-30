function emdlab_flib_setPlot()

fig = gcf;
ax = gca;

ax.XLabel.FontName = 'Helvetica';
ax.XLabel.FontSize = 14;
ax.XAxis.FontSize = 12;

ax.YLabel.FontName = 'Helvetica';
ax.YLabel.FontSize = 14;
ax.YAxis.FontSize = 12;

fig.Position = [0,0,600,300];
movegui(fig, 'center');
end