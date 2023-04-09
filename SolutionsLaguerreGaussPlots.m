
addpath("PlotPub\lib")
clear
close all

get_default_figure()


x = linspace(0,30,1024);

n_max = 7;
m_max = 7;

for n_value = 1:n_max
    for m_value = 1:m_max
        if m_value <=n_value

            Ln  = LaguerreG(n_value,m_value,x);
            XLn = XLaguerreG(n_value,m_value,x);
            H = Ln+1i*XLn;
            fig = figure(1);
%             fig.Position = [2429 367 829 399];
            plot(x,Ln)
            hold on
            plot(x,XLn)
            plot(x,abs(H))
            yline(0,'-.');
            hold off
            lg1 = sprintf('L_{%d}^{%d}(x)', n_value,m_value);
            lg2 = sprintf('X_{%d}^{%d}(x)', n_value,m_value);
            lg3 = sprintf('HL_{%d}^{%d}(x)', n_value,m_value);
            legend(lg1,lg2,lg3)
            xlim([0,3*(n_value+1)+m_value])
            ylim([-3*(m_value*n_value)-2*n_value*(2^n_value),3*(m_value*n_value)+2*n_value*(2^n_value)])
            set(gca, 'FontSize', 14)
            path_file = ['LaguerreGSol_n',num2str(n_value),'_m',num2str(m_value),'.png'];
            xlabel('x')
            print(path_file, '-dpng', '-r600')
        end
    end
end





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
