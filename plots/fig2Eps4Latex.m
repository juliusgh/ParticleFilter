%% Matlab-Script: fig2Eps4Latex.m  -  Set labels to latex, by J. Stoerkle

%# create a test plot:
% figure;plot(0:0.1:10,sin(0:0.1:10),'linewidth',1.5,'displayname','sin$(\Omega t)$'); xlabel('a test $x$');ylabel('test $y$'); legend('show')

%# set latex to the most important labels
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.XLabel.Interpreter   = 'latex';
ax.XLabel.FontSize      = 12;
ax.YLabel.Interpreter   = 'latex';
ax.YLabel.FontSize      = 12;
if ~isempty(ax.Legend)
    ax.Legend.Interpreter = 'latex';
    ax.Legend.FontSize    = 12;
end
ax.Title.Interpreter = 'latex';
ax.FontSize = 12;
box on
grid on

%% Plot results
set(gcf,'PaperPositionMode','auto') % take print-pixel from window-pixel
set(gcf,'units','pixels')           % also possibel: 'centimeters'
pos = get(gcf,'Position');

%# set custom heigth and width
% myWidthPx = 570; % Diss:  \includegraphics[scale=1.0]{Bilder/myFig.eps}
myWidthPx = 517; % Studi: \includegraphics[scale=1.0]{Bilder/myFig.eps}
set(gcf,'Position',[pos(1) pos(2) myWidthPx*0.5 myWidthPx*0.5]); % 450 270

%# set minimal white space 
%# (https://de.mathworks.com/help/matlab/creating_plots/save-figure-with-minimal-white-space.html)
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

%% save as eps
print(gcf,'-depsc','-loose',['-r','200'],'myFig'); % since R2014