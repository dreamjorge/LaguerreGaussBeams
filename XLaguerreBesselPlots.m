
addpath("PlotPub\lib")
clear
% close all

get_default_figure()

n_points = 2^10;
x = linspace(0,10,n_points);

% n_value = 50;
% m_value = 11;
% n_value = 13;
% m_value = 2;

LGnm = XLaguerreG(n_value,m_value,x);
% bessely(nu,Z)
N = n_value + (m_value+1)/2;
J = 1.... (factorial(n_value+m_value)/factorial(m_value))...
    ...*exp(x/2)...
    ...*(N.*x).^(-m_value/2)...
    .*bessely(m_value,2*sqrt(N*x));
fig = figure(1);
%             fig.Position = [2429 367 829 399];
plot(x,LGnm)
hold on
plot(x,J,"--")
yline(0,'-.');
hold off
lg1 = '$B_{mn} {\textbf{e}}^{-x/2}x^{m/2}X_{n}^{m}(x)$';
lg2 = '$N_m(2\sqrt{Nx})$';
legend(lg1,lg2,'Interpreter','latex')
...,'interpreter','latex')
set(gca, 'FontSize', 14)
path_file = ['XLaguerreNBessel_n',num2str(n_value),'_m',num2str(m_value),'.png'];
xlabel('x')
xlim([0,10])
ylim([-0.6, 0.6])
print(path_file, '-dpng', '-r600')







function get_default_figure()

% Defaults for this blog post
width = 8;     % Width in inches
height = 5;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 30;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize

% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
% set(0,'FontSize',     fsz);
% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1)-40, defpos(2)-40, width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

end
