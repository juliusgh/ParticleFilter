function itm_formatfig(varargin)
% itm_formatfig(varargin)
%
% The script itm_formatfig formats Matlab figures, it is suitable for
% standard xy-plots and for figures with subplots. For other figures such
% as 3D plots or bar diagrams, a basic formation is carried out.
% If desired, figures are stored as fig-, pdf- and png-documents.
%
% By default, this script will apply to the currently active figure.
% It is advisable to carry out the formatting directly after your figure
% has been created, e.g.
%
% figure;
% plot([0:1e-3:2], sin(2*pi*[0:1e-3:2]));
% xlabel('time [s]');
% ylabel('displacement [m]');
% grid on;
% itm_formatfig(2);
%
%
% SELECT THE PRESET FOR THE DESIRED FIGURE FORMAT
%
% The first input argument is obligatory and defines the figure format preset as follows.
%
%
% 1 or 'LatexNarrow'..... Creates a narrow figure for Latex documents,
%                         fits e.g. the ITM Latex template for BSC, MSC, ...
%                         Figure size is 6 cm x 8 cm.
% 2 or 'LatexWide' ...... Creates a wide figure for Latex documents,
%                         fits e.g. the ITM Latex template for BSC, MSC, ...
%                         Figure size is 14 cm x 8 cm.
%                         This is the default format.
% 3 or 'Word' ........... Creates a figure for Word documents
%                         Figure size is 13 cm x 10 cm.
% 4 or 'PowerPoint' ..... Creates a figure for Powerpoint presentations,
%                         please use .png format!
% 5 or 'Poster' ......... Sets the aspect ratio to 4 x 3, here the size is
%                         12 cm x 9 cm
% 6 or 'WikiNarrow' ..... Creates a narrow figure for the Wiki
% 7 or 'WikiWide' ....... Creates a wide figure for the Wiki
% 8 or 'Movie' .......... Sets the aspect ratio to 4 x 3, here the size is
%                         12 cm x 9 cm
% 9 or 'UserDef' ........ User-definded preset
% 10 or 'Diss' .......... Creates a figure that fits the ITM Latex template
%                         for dissertations
%
% Example:
% itm_formatfig(2);
%
%
% STORE THE FORMATTED FIGURE AND EXPORT TO A GRAPHICS FORMAT
%
% In order to store the formatted figure, the parameter/value pairs
%
% 'FigName', '<filename>' for the filename and optionally
% 'FigPath', '<directory>' for the target directory.
%
% can be passed. Files of the types .fig, .png and .pdf will be created and stored.
% In case no target directory is defined, the current Matlab working directory will be used.
%
% Example:
% itm_formatfig(2,'FigName','mypicture','FigPath','~');
%
%
% FURTHER FIGURE PARAMETERS
%
% Apart from the presets, other input options may be defined as parameter/value pairs as follows
% (default values are given in {}).
%
% 'BorderWidth' ......... Width of box border lines [pt] (numeric scalar)
% 'Box' ................. Display a box or not ({'on'}, 'off)
% 'FigName'...............Filename for saving fig, pdf and png documents in
%                         the current working directory (string)
% 'FigPath'...............Target directory if other than the current working
%                         directory is desired (string).
% 'FigSize' ............. Figure width and height [cm] (numeric 2-by-1 vector)
% 'Figure' .............. Figure can be specified ({gcf}, figure handle)
% 'FontName' ............ Font to be used in all contained text, e.g. 'Arial'
%                         (string)
% 'FontSize' ............ Fontsize in all text [pt] (numeric scalar)
% 'GridLineStyle' ....... Style for the grid lines ('-', '--', {:}, '-.')
% 'Grid' ................ Switch the grid on or off ({'on'}, 'off')
% 'LabelFlag' ........... set 0 to remove label text and title from the
%                         figure ({1}, 0)
% 'LegendLocation' ...... Location of the legend ({'NorthEast'}, 'NorthWest',
%                         'SouthEastOutside' etc.)
% 'LineWidth' ........... Linewidth for all lines [pt] (numeric scalar)
% 'MarkerSize' .......... Markersize for all lines. If empty, all marker
%                         sizes are kept (numeric scalar)
% 'SpreadSubplots' ...... Larger subplot diagrams within the given figure
%                         size. Specify the number of subplots as a numerical 2-by-1
%                         vector in vertical and horizontal direction. {[]}
% 'OnlyBasic' ........... Logical parameter, whether to quit after the very
%                         basic formatting. This avoids problems with 3D
%                         plots and subplots. When in use, the adjustment
%                         of the position of the axes in the figure as well
%                         as the positioning of the labels is omitted.
%                         ('true', {'false'})
% 'Verbose' ..............Switch text output to the matlab console on or off.
%                         ({'on'}, 'off')
%
%
% Examples:
% itm_formatfig(2,'LegendLocation','NorthWestOutside','FigName','mypicture');
% itm_formatfig('Word','LineWidth',2,'FontName','Times','LegendLocation','SouthEast');
% itm_formatfig(2,'SpreadSubplots', [2 2]);
% itm_formatfig(2,'LineWidth',2,'FontName','Arial','FigName','bild1','figpath','~');
% itm_formatfig('PowerPoint','LegendLocation','NorthEast','FigName','mypicture','FigPath','~');
%
%

% get matlab version
ver_ = version('-release');
if str2double(ver_(1:4))>2013
    newMatlab_ = 1;
else
    newMatlab_ = 0;
end

% Select default settings
fontsize_ = 12;
fontname_ = 'Times';
linewidth_ = 0.5;
markerSize_ = 5;
epsxsize_   = 13.5;
epsysize_   = 8;
fixedWhiteSpace_ = false;
usebox_     = 'on';
boundaryWidth_ = [1 1 1 1];
borderwidth_ = 0.5;
gridlinestyle_ = ':';
switchgrid_ = 'on';
legendlocation_ = 'North';
labelflag_ = 1;
verbose_ = 1; % parameter to adjust output, 0: off
backcolor_ = [1 1 1]; % white background
legendbox_ = 'on'; % box around the legend
spreadSubplots_ = [];
onlybasic_ = false; % Perform only the very basic adjustments, supposed to work always
figHandle_ = gcf;
setDefaultAxesColorOrder_ = false;
figname_ = 0;
figpath_ = './';

% Treat input arguments
if(nargin == 0)
    error('itm_formatfig: Please define format preset as first input argument, e.g. itm_formatfig(''LatexWide''). The script itm_formatfig is aborted.')
    return;
end;


% Check verbose
h_ = [];
if length(varargin) > 1
    for h_ = 2:2:length(varargin)
        if strcmp(lower(varargin{h_}),'verbose')
            
            if ~isa(varargin{h_+1},'char')
                error('itm_formatfig: Value of ''Verbose'' must be either ''on'' or ''off''. The script itm_formatfig is aborted.');
                return;
            else
                if ~strcmp(varargin{h_+1},'on') && ~strcmp(varargin{h_+1},'off')
                    error('itm_formatfig: Value of ''Verbose'' must be either ''on'' or ''off''. The script itm_formatfig is aborted.');
                    return;
                end
            end
            
            if strcmp(varargin{h_+1},'on')
                verbose_ = 1;
            end
            if strcmp(varargin{h_+1},'off')
                verbose_ = 0;
            end
            
        end
    end
end
h_ = [];



% Presets
switch lower(varargin{1})
    case {1,'1', 'latexnarrow'}
        epsxsize_   = 6;
        epsysize_   = 6;
        linewidth_  = 0.5;
        borderwidth_ = 0.5;
        if(verbose_>0)
            fprintf('itm_formatfig: Preset 1 ''LatexNarrow''\n');
        end
    case {2,'2','latexwide'}
        if(verbose_>0)
            fprintf('itm_formatfig: Preset 2 ''LatexWide''\n');
        end
    case {3,'3','word'}
        epsxsize_ = 13;
        epsysize_ = 10;
        fontname_ = 'Arial';
        fontsize_ = 11;
        linewidth_ = 1;
        gridlinestyle_ = ':';
        borderwidth_ = 1;
        legendlocation_ = 'NorthEast';
        if(verbose_>0)
            fprintf('itm_formatfig: Preset 3 ''Word''\n');
        end
    case {4,'4','pp','powerpoint'}
        epsxsize_ = 11.5;
        epsysize_ = 8;
        fontsize_ = 16;
        fontname_ = 'Arial';
        linewidth_ = 1.5;
        borderwidth_ = 1;
        markerSize_ = 6;
        if(verbose_>0)
            fprintf('itm_formatfig: Preset 4 ''PowerPoint''\n');
        end
    case {5,'5','poster'}
        epsxsize_ = 40;
        epsysize_ = 20;
        fontsize_ = 21;
        linewidth_ = 2;
        markerSize_ = 8;
        borderwidth_ = 1;
        if(verbose_>0)
            fprintf('itm_formatfig: Preset 5 ''Poster''\n');
        end
    case {6,'6','wikinarrow'}
        epsxsize_   = 8;
        epsysize_   = 6;
        fontsize_ = 10;
        linewidth_ = 1;
        if(verbose_>0)
            fprintf('itm_formatfig: Preset 6 ''WikiNarrow''\n');
        end
    case {7,'7','wikiwide'}
        epsxsize_   = 12;
        epsysize_   = 6;
        fontsize_ = 10;
        linewidth_ = 1;
        if(verbose_>0)
            fprintf('itm_formatfig: Preset 7 ''WikiWide''\n');
        end
    case {8,'8','movie'}
        epsxsize_   = 12;
        epsysize_   = 9;
        fontsize_ = 10;
        if(verbose_>0)
            fprintf('itm_formatfig: Preset 8 ''Movie''\n');
        end
    case {9,'9','userdef'}
        if(verbose_>0)
            fprintf('itm_formatfig: Preset 9 ''UserDef''\n');
        end
    case {10,'10','diss'}
        epsxsize_   = 15.99;
        epsysize_   = 9.88;
        fontsize_ = 11;
        if(verbose_>0)
            fprintf('itm_formatfig: Preset 10 ''Diss''\n');
        end
        
    otherwise
        %             fprintf('\nFormat preset (first input argument) not defined or invalid. Default format "2: latexwide" is used.\n');
        error('itm_formatfig: Format preset is not defined or is invalid, please check first input argument. The script itm_formatfig is aborted.')
        return;
end


if length(varargin) > 1 && rem(length(varargin),2) == 0
    error('itm_formatfig: Number of input arguments is invalid, please check parameter/value pairs. The script itm_formatfig is aborted.')
    return;
end

idx_ = 3;

% Treat optional input arguments
for h_ = idx_:2:length(varargin)
    
    switch lower(varargin{h_-1})
        case 'axes'
            if(~isempty(varargin{h_}))
                axes(varargin{h_});
            end
        case 'borderwidth'
            if(~isempty(varargin{h_}))
                if ~isa(varargin{h_},'numeric')
                    error('itm_formatfig: Value of ''BorderWidth'' must be a numeric scalar. The script itm_formatfig is aborted.');
                    return;
                else
                    size_ = [];
                    size_ = size(varargin{h_});
                    if size_(1) > 1 || size_(2) >1
                        error('itm_formatfig: Value of ''BorderWidth'' must be a numeric scalar. The script itm_formatfig is aborted.');
                        clear size_;
                        return;
                    end
                end
                borderwidth_ = varargin{h_};
            end
        case 'boundarywidth'
            if(~isempty(varargin{h_}))
                boundaryWidth_ = varargin{h_};
                fixedWhiteSpace_ = true;
            end
        case 'box'
            if(~isempty(varargin{h_}))
                
                if ~isa(varargin{h_},'char')
                    error('itm_formatfig: Value of ''Box'' must be either ''on'' or ''off''. The script itm_formatfig is aborted.');
                    return;
                else
                    if ~strcmp(varargin{h_},'on') && ~strcmp(varargin{h_},'off')
                        error('itm_formatfig: Value of ''Box'' must be either ''on'' or ''off''. The script itm_formatfig is aborted.');
                        return;
                    end
                end
                
                usebox_ = varargin{h_};
            end
        case 'defaultaxescolororder'
            if(~isempty(varargin{h_}))
                setDefaultAxesColorOrder_ = logical(varargin{h_});
            end
        case 'figsize'
            if(~isempty(varargin{h_}))
                if ~isa(varargin{h_},'numeric')
                    error('itm_formatfig: Value of ''FigSize'' must be a numeric 2-by-1 vector. The script itm_formatfig is aborted.');
                    return;
                else
                    if size(varargin{h_}, 1) > 2 || size(varargin{h_}, 2) ~= 1
                        error('itm_formatfig: Value of ''FigSize'' must be a numeric 2-by-1 vector. The script itm_formatfig is aborted.');
                        clear size_;
                        return;
                    end
                end
                
                epsxsize_ = varargin{h_}(1);
                epsysize_ = varargin{h_}(2);
            end
        case 'figure'
            if(~isempty(varargin{h_}))
                figHandle_ = varargin{h_};
                figure(figHandle_(1));
            end
        case 'fixedwhitespace'
            if(~isempty(varargin{h_}))
            elseif(islogical(varargin{h_}))
                fixedWhiteSpace_ = varargin{h_};
            elseif(isnumeric(varargin{h_}))
                fixedWhiteSpace_ = true;
                boundaryWidth_ = varargin{h_};
            elseif(ischar(varargin{h_}))
                if(strcmpi('true',varargin{h_}))
                    fixedWhiteSpace_ = true;
                elseif(strcmpi('false',varargin{h_}))
                    fixedWhiteSpace_ = false;
                else
                    error('Unknown argument ''%s'' passed together with ''%s''!',varargin{h_},varargin{h_-1});
                end
            else
                error('Unknown argument ''%s'' passed together with ''%s''!',varargin{h_},varargin{h_-1});
            end
        case 'fontname'
            if(~isempty(varargin{h_}))
                
                
                if ~isa(varargin{h_},'char')
                    error('itm_formatfig: Value of ''FontName'' must be a string. The script itm_formatfig is aborted.');
                    return;
                end
                
                fontname_ = varargin{h_};
            end
        case 'fontsize'
            if(~isempty(varargin{h_}))
                if ~isa(varargin{h_},'numeric')
                    error('itm_formatfig: Value of ''FontSize'' must be a numeric scalar. The script itm_formatfig is aborted.');
                    return;
                else
                    if size(varargin{h_}, 1) > 1 || size(varargin{h_}, 2) > 1
                        error('itm_formatfig: Value of ''FontSize'' must be a numeric scalar. The script itm_formatfig is aborted.');
                        return;
                    end
                end
                fontsize_ = varargin{h_};
            end
        case 'grid'
            if(~isempty(varargin{h_}))
                if ~isa(varargin{h_},'char')
                    error('itm_formatfig: Value of ''Grid'' must be either ''on'' or ''off''. The script itm_formatfig is aborted.');
                    return;
                else
                    if ~strcmp(varargin{h_},'on') && ~strcmp(varargin{h_},'off')
                        error('itm_formatfig: Value of ''Grid'' must be either ''on'' or ''off''. The script itm_formatfig is aborted.');
                        return;
                    end
                end
                switchgrid_ = varargin{h_};
            end
        case 'gridlinestyle'
            if(~isempty(varargin{h_}))
                gridlinestyle_ = varargin{h_};
            end
        case 'labelflag'
            if(~isempty(varargin{h_}))
                labelflag_ = varargin{h_};
            end
        case 'legendbox'
            if(~isempty(varargin{h_}))
                
                if ~isa(varargin{h_},'char')
                    error('itm_formatfig: Value of ''LegendBox'' must be either ''on'' or ''off''. The script itm_formatfig is aborted.');
                    return;
                else
                    if ~strcmp(varargin{h_},'on') && ~strcmp(varargin{h_},'off')
                        error('itm_formatfig: Value of ''LegendBox'' must be either ''on'' or ''off''. The script itm_formatfig is aborted.');
                        return;
                    end
                end
                legendbox_ = varargin{h_};
            end
        case 'legendlocation'
            if(~isempty(varargin{h_}))
                legendlocation_ = varargin{h_};
            end
        case 'linewidth'
            % Empty keeps all corresponding settings
            if ~isa(varargin{h_},'numeric')
                error('itm_formatfig: Value of ''LineWidth'' must be a numeric scalar. The script itm_formatfig is aborted.');
                return;
            else
                size_ = [];
                size_ = size(varargin{h_});
                if size_(1) > 1 || size_(2) >1
                    error('itm_formatfig: Value of ''LineWidth'' must be a numeric scalar. The script itm_formatfig is aborted.');
                    clear size_;
                    return;
                end
            end
            linewidth_ = varargin{h_};
        case 'markersize'
            % Empty keeps all corresponding settings
            
            if ~isa(varargin{h_},'numeric')
                error('itm_formatfig: Value of ''MarkerSize'' must be a numeric scalar. The script itm_formatfig is aborted.');
                return;
            else
                size_ = [];
                size_ = size(varargin{h_});
                if size_(1) > 1 || size_(2) >1
                    error('itm_formatfig: Value of ''MarkerSize'' must be a numeric scalar. The script itm_formatfig is aborted.');
                    clear size_;
                    return;
                end
            end
            
            markerSize_ = varargin{h_};
        case 'spreadsubplots'
            % Empty keeps all corresponding settings
            spreadSubplots_ = varargin{h_};
        case 'onlybasic'
            if(~isempty(varargin{h_}))
                onlybasic_ = varargin{h_};
            end
            
            
        case 'figname'
            if(~isempty(varargin{h_}))
                
                if ~isa(varargin{h_},'char')
                    error('itm_formatfig: Value of ''FigName'' must be a string. The script itm_formatfig is aborted.');
                    return;
                end
                
                figname_ = varargin{h_};
            end
        case 'figpath'
            if(~isempty(varargin{h_}))
                
                
                if ~isa(varargin{h_},'char')
                    error('itm_formatfig: Value of ''FigPath'' must be a string. The script itm_formatfig is aborted.');
                    return;
                end
                
                figpath_ = varargin{h_};
            end
            
        case 'verbose'
            % was checked above
            
            
        otherwise
            error('itm_formatfig: Unknown option ''%s'' is passed. The script itm_formatfig is aborted.',varargin{h_-1});
            
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set defaults
% set(0,'DefaultAxesFontsize'   ,fontsize_)
% set(0,'DefaultTextFontsize'   ,fontsize_)
set(0,'DefaultTextInterpreter','tex');
set(0,'DefaultTextLineStyle'  ,'none');
if(setDefaultAxesColorOrder_)
    set(0,'DefaultAxesColorOrder',[0 0 1;1 0 0;0 0 0;0 0.5 0],...
        'DefaultAxesLineStyleOrder','-|-|-|-');
end

% store Renderer
set(gcf,'Renderer','painters');
% renderer_ = get(gcf,'Renderer');


% Ensure the windowstyle is not 'docked' as then some things cannot be
% changed
set(gcf,'WindowStyle','normal');

mainaxeshandles_   = [];
legendaxeshandles_ = [];
% get handles
linehandles_ = findobj(gcf,'Type','Line');
if newMatlab_
    allaxeshandles_ = findobj(gcf,'type','axes');
    mainaxeshandles_ = findobj(gcf,'type','axes');
    if ~isempty(findobj(gcf,'type','legend'))
        legendaxeshandles_ = findobj(gcf,'type','legend');
    end
    
else
    allaxeshandles_ = findobj(gcf,'type','axes');
    % Find all axes and legends
    handletype = strcmp('', get(allaxeshandles_, 'Tag'));
    for h = 1 : length(allaxeshandles_)
        if handletype(h)
            mainaxeshandles_(end+1) = allaxeshandles_(h);
        else
            legendaxeshandles_(end+1) =  allaxeshandles_(h);
        end
    end
end

% If a linewidth is set, adjust the width of all lines
if(~isempty(linewidth_))
    set(linehandles_,'LineWidth',linewidth_);
end
if(~isempty(markerSize_))
    set(linehandles_,'MarkerSize',markerSize_);
end
% Adjust axes
% switch box on/off
set(mainaxeshandles_, 'Box', usebox_);
set(legendaxeshandles_, 'Box', legendbox_);


% adjust borderwidth of box
if newMatlab_
    set(mainaxeshandles_,'LineWidth', borderwidth_);
    set(legendaxeshandles_,'LineWidth', borderwidth_);
else
    set(allaxeshandles_,'LineWidth', borderwidth_);
end

% set grid
set(mainaxeshandles_, 'GridLineStyle', gridlinestyle_);
for h = 1:length(mainaxeshandles_)
    grid(mainaxeshandles_(h), switchgrid_);
end

% format text
alltext_ = [findobj(gcf,'Type','text'); findobj(gcf,'Tag','legend')];
myaxes_ = findobj(gcf,'Type','axes');
for h_ = 1:length(myaxes_)
    alltext_ = [alltext_; myaxes_(h_); get(myaxes_(h_),'Title'); ...
        get(myaxes_(h_),'XLabel'); get(myaxes_(h_),'YLabel'); ...
        get(myaxes_(h_),'ZLabel')];
end
set(alltext_,'FontSize',fontsize_);
set(alltext_,'FontName',fontname_);
set(gca,'FontName', fontname_);

% get fontheight
fontheight = fontsize_ * 127/3600; % fontheight in cm


% remove all text, if desired
if labelflag_ == 0;
    set(gca,'XTickLabel',[],'XTickLabelMode','manual');% remove x ticks
    set(gca,'YTickLabel',[],'YTickLabelMode','manual');% remove y ticks
    title('');
    set(gca,'Position',[.01 .01 .98 .98]);             % 98 % of length and height
end

% prepare figure
set(gcf,'PaperType','a4');                             % paper size A4
set(gcf,'PaperOrientation','portrait');                % orientation portrait
set(gcf,'PaperUnits','centimeters');                   % paper unit Centimeters
set(gcf,'PaperPosition', [4 4 epsxsize_ epsysize_]);   % size of figure
set(gcf,'PaperPositionMode','manual');
set(gcf,'Units','centimeters');                        % figure unit Centimeters
fp = get(gcf,'Position');                              % get figure position
set(gcf,'Position',[fp(1) fp(2) epsxsize_ epsysize_]); % make figure wysiwyg

% Farben richten. Dabei wird vorausgesetzt, dass die Axis durchscheinend
% ist
% Dies sollte im startup durchgefuehrt werden.
set (gcf,'Color',backcolor_);   % set colour of figure
% set (gca,'Color',backcolor);    % set colour of axis
set(gcf,'InvertHardCopy','off') % Damit's auch auf dem Drucker so rauskommt
% Adjust legend location
if newMatlab_
    set(findobj(gcf,'Type','legend'),'Location', legendlocation_);
else
    set(findobj(gcf,'Type','axes','Tag','legend'),'Location', legendlocation_);
end

% Determine whether the selected figure is suitable for full formatting or
% if a basic version shall be used, which does not break anything
[az_,el_] = view;

% Check if there is a colorbar

Colorbar_flag = 0;
ax = findall(gcf,'type','axes');
axtags = get(ax,'Tag');
Colorbar_flag = ismember('Colorbar', axtags);

% ab MatlabR2014b
if Colorbar_flag == 0
    ax_neu = findall(gcf,'type','ColorBar');
    Colorbar_flag = ~isempty(ax_neu);
end

if(length(mainaxeshandles_) > 1 && isempty(spreadSubplots_))
    [az_,el_] = view;
    % Subplots in use, do only basic formatting
    onlybasic_ = true;
    %fprintf('\nFound subplots, quitting after basic formatting\n')
    fprintf('Found subplots, basic formatting.\n')
elseif(az_~=0 || el_~=90)
    % Not a 2D plot but some 3D graph
    onlybasic_ = true;
    %fprintf('\nFound 3D-Plot, quitting after basic formatting\n')
    fprintf('Found 3D-Plot, basic formatting.\n')
elseif(Colorbar_flag == 1)
    % Not a 2D plot but some 3D graph
    onlybasic_ = true;
    % fprintf('\nFound Colorbar, quitting after basic formatting\n')
    fprintf('Found Colorbar, basic formatting.\n')
end

if(onlybasic_)
    % Don't adjust sizes or positions, leave everything to user
    % Set units of all axes to normalized, this enables the manual resizing
    if newMatlab_
        set(mainaxeshandles_,'Units','normalized');
        set(legendaxeshandles_,'Units','normalized');
    else
        set(allaxeshandles_,'Units','normalized');
    end
    
    
    % set(allaxeshandles_,'Units','normalized');
    
    %     set(gcf,'renderer',renderer_); % Reset renderer
    % Usage hints
    
    % optional saving
    if figname_ ~= 0
        if figpath_(end) ~= '/'
            figpath_ = [figpath_ '/'];
        end
        saveas(gcf,[figpath_ figname_ '.fig']);
        print(gcf,'-depsc',[figpath_ figname_ '.eps']);
        print(gcf,'-dpng',[figpath_ figname_ '.png']);
    else
        if(verbose_>0)
            fprintf('\nTo store the figure, please use one of these commands');
            fprintf('\n\tprint -depsc ''myFigure.eps''\n\tprint -dpng  ''myFigure.png''\n');
        end
    end
    
    
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% format figure
if(length(mainaxeshandles_) == 1)
    %% One axes, general adjustments
    % If only one axes is available, adjust settings
    xlabelh = get(mainaxeshandles_, 'XLabel');
    ylabelh = get(mainaxeshandles_, 'YLabel');
    % Set units
    set(xlabelh,'Units','centimeters');
    set(ylabelh,'Units','centimeters');
    set(mainaxeshandles_,'Units','centimeters');
    %  axpos = get(mainaxeshandles_,'Position');
    % initialize label position
    set(xlabelh,'Position',[0 0]);
    set(ylabelh,'Position',[0 0]);
    % get tight inset of main axes
    ti = get(mainaxeshandles_, 'TightInset');
    
    % check if ticklabels were overwritten by text
    flag_x = zeros(1,length(alltext_));
    flag_y = zeros(1,length(alltext_));
    for k_ = 1 : length(alltext_)
        flag_x(k_) = strcmp('MUXTL', get(alltext_(k_), 'Tag'));
        flag_y(k_) = strcmp('MUYTL', get(alltext_(k_), 'Tag'));
    end
    if ((norm(flag_x) > 0) || (norm(flag_y) > 1))
        ti(1) = ti(1) + 1.6*fontheight;
        ti(2) = ti(2) + 1.2*fontheight;
    end
    
    % initialize new label positions
    posxlabel = get(xlabelh,'Position');
    posylabel = get(ylabelh,'Position');
    
    % Adjust xlabel position
    set(xlabelh,'Margin',1);
    set(xlabelh,'VerticalAlignment','bottom');
    % set  vertical position
    %     posxlabel(2) = -ti(2) -1.5*fontheight;
    posxlabel(2) = -ti(2) -1.2*fontheight;
    set(xlabelh,'Position', posxlabel);
    
    % Adjust ylabel position
    set(ylabelh,'Margin',1);
    set(ylabelh,'VerticalAlignment','top');
    % set horizontal position
    posylabel(1) = - ti(1) -1.6*fontheight;
    set(ylabelh,'Position',posylabel);
    
    % Adjust position of axes with respect to the labels
    axpos = get(mainaxeshandles_,'position');
    posxlabel = get(xlabelh,'Position');
    posylabel = get(ylabelh,'Position');
    
    axpos(1) = -posylabel(1) + 0.35*fontheight;
    axpos(2) = -posxlabel(2) + 0.3*fontheight;
    %         axpos(1) = -posylabel(1);
    %         axpos(2) = -posxlabel(2);
    
    
    set(mainaxeshandles_,'position',axpos);
    
    %% One axes, legend specific adjustments
    if ~isempty(legendaxeshandles_)
        % Legends in the figure
        if strfind(legendlocation_, 'Outside') > 0
            % legend outside
            set(gcf,'Units','centimeters');
            set(mainaxeshandles_,'Units','centimeters');
            set(legendaxeshandles_,'Units','centimeters');
            
            % Move legend to the boundary, adjust the size of the axes so it fills the space
            figpos = get(gcf,'position');
            legendpos = get(legendaxeshandles_,'position');
            legendpos(1) = figpos(3) - legendpos(3);
            set(legendaxeshandles_,'position',legendpos);
            axpos = get(mainaxeshandles_,'position');
            axpos(3) = figpos(3) - axpos(1) - legendpos(3) - ti(3) - 0.5*fontheight;
            set(mainaxeshandles_,'Position',axpos);
            axpos = get(mainaxeshandles_,'Position');
            axpos(4) = figpos(4) - axpos(2) - ti(4);
            set(mainaxeshandles_, 'Position',axpos);
            
            axpos = get(mainaxeshandles_,'position');
            posxlabel = get(xlabelh,'Position');
            posylabel = get(ylabelh,'Position');
            
            axpos(1) = -posylabel(1) + 0.35*fontheight;
            axpos(2) = -posxlabel(2) + 0.3*fontheight;
            set(mainaxeshandles_,'position',axpos);
            
            % Adjust position of the legend
            legendpos = get(legendaxeshandles_,'Position');
            legendpos(2) = axpos(2) + axpos(4) - legendpos(4);
            set(legendaxeshandles_,'Position',legendpos);
            
            axpos = get(mainaxeshandles_,'position');
            axpos(3) = figpos(3) - axpos(1) - legendpos(3) - ti(3);
            set(mainaxeshandles_,'Position',axpos);
            
            axpos = get(mainaxeshandles_,'position');
            posxlabel = get(xlabelh,'Position');
            posylabel = get(ylabelh,'Position');
            
            axpos(1) = -posylabel(1)+ 0.3*fontheight;
            axpos(2) = -posxlabel(2)+ 0.3*fontheight;
            set(mainaxeshandles_,'position',axpos);
            
        else
            % legend inside
            set(gcf,'Units','centimeters');
            set(mainaxeshandles_,'Units','centimeters');
            set(legendaxeshandles_,'Units','centimeters');
            figpos = get(gcf,'position');
            axpos = get(mainaxeshandles_,'position');
            axpos(3) = figpos(3) - axpos(1) - ti(3) - 0.5*fontheight;
            set(mainaxeshandles_,'Position',axpos);
            axpos = get(mainaxeshandles_,'Position');
            axpos(4) = figpos(4) - axpos(2) - ti(4);
            set(mainaxeshandles_, 'Position',axpos);
            
            axpos = get(mainaxeshandles_,'position');
            posxlabel = get(xlabelh,'Position');
            posylabel = get(ylabelh,'Position');
            
            axpos(1) = -posylabel(1)  + 0.3*fontheight;
            axpos(2) = -posxlabel(2)  + 0.3*fontheight;
            set(mainaxeshandles_,'position',axpos);
            
            set(legendaxeshandles_, 'Location', legendlocation_);
            
            axpos = get(mainaxeshandles_,'position');
            posxlabel = get(xlabelh,'Position');
            posylabel = get(ylabelh,'Position');
            
            axpos(1) = -posylabel(1)  + 0.3*fontheight;
            axpos(2) = -posxlabel(2)  + 0.3*fontheight;
            set(mainaxeshandles_,'position',axpos);
            
        end
        % Center labels, change and reset all units
        set(gcf,'Units','normalized');
        set(mainaxeshandles_,'Units','normalized');
        set(xlabelh,'Units','normalized');
        set(ylabelh,'Units','normalized');
        
        set(xlabelh,'HorizontalAlignment','center');
        posxlabel = get(xlabelh,'Position');
        posxlabel(1) = 0.5;
        set(xlabelh,'Position',posxlabel);
        
        set(ylabelh,'HorizontalAlignment','center');
        posylabel = get(ylabelh,'Position');
        posylabel(2) = 0.5;
        set(ylabelh,'Position',posylabel);
        
        set(gcf,'Units','centimeters');
        set(mainaxeshandles_,'Units','centimeters');
        set(xlabelh,'Units','centimeters');
        set(ylabelh,'Units','centimeters');
        
        
        % If the white space arouns the figure is to be fixed, handle that
        if(fixedWhiteSpace_)
            
            oldPos_ = get(mainaxeshandles_,'Position');
            
            newPos_(3) = oldPos_(3) - boundaryWidth_(1) - boundaryWidth_(3);
            newPos_(4) = oldPos_(4) - boundaryWidth_(2) - boundaryWidth_(4);
            newPos_(1:2) = oldPos_(1:2) + boundaryWidth_(1:2);
            
            set(mainaxeshandles_,'Position',newPos_);
            
            set(gcf,'Units','normalized');
            set(mainaxeshandles_,'Units','normalized');
            set(xlabelh,'Units','normalized');
            set(ylabelh,'Units','normalized');
            
            set(xlabelh,'HorizontalAlignment','center');
            posxlabel = get(xlabelh,'Position');
            posxlabel(1) = 0.5;
            set(xlabelh,'Position',posxlabel);
            
            set(ylabelh,'HorizontalAlignment','center');
            posylabel = get(ylabelh,'Position');
            posylabel(2) = 0.5;
            set(ylabelh,'Position',posylabel);
            
            set(gcf,'Units','centimeters');
            set(mainaxeshandles_,'Units','centimeters');
            set(xlabelh,'Units','centimeters');
            set(ylabelh,'Units','centimeters');
            
            if strfind(legendlocation_, 'Outside') > 0
                legendpos(1) = legendpos(1) - boundaryWidth_(3);
                legendpos(2) = legendpos(2) - boundaryWidth_(4);
                set(legendaxeshandles_,'Position',legendpos);
            end
        end
        
    else
        %% One axes, no legend
        set(gcf,'Units','centimeters');
        set(mainaxeshandles_,'Units','centimeters');
        figpos = get(gcf,'position');
        axpos = get(mainaxeshandles_,'position');
        axpos(3) = figpos(3) - axpos(1) - ti(3) - 0.5*fontheight;
        set(mainaxeshandles_,'Position',axpos);
        axpos = get(mainaxeshandles_,'Position');
        axpos(4) = figpos(4) - axpos(2) - ti(4);
        set(mainaxeshandles_, 'Position',axpos);
        
        %         axpos = get(mainaxeshandles_,'position');
        %         posxlabel = get(xlabelh,'Position');
        %         posylabel = get(ylabelh,'Position');
        %
        %         axpos(1) = -posylabel(1);
        %         axpos(2) = -posxlabel(2);
        %         set(mainaxeshandles_,'position',axpos);
        
        axpos = get(mainaxeshandles_,'position');
        posxlabel = get(xlabelh,'Position');
        posylabel = get(ylabelh,'Position');
        
        axpos(1) = -posylabel(1) + 0.35*fontheight;
        axpos(2) = -posxlabel(2) + 0.3*fontheight;
        %         axpos(1) = -posylabel(1);
        %         axpos(2) = -posxlabel(2);
        set(mainaxeshandles_,'position',axpos);
        
    end
    % Center labels a second time
    set(gcf,'Units','normalized');
    set(mainaxeshandles_,'Units','normalized');
    set(xlabelh,'Units','normalized');
    set(ylabelh,'Units','normalized');
    
    set(xlabelh,'HorizontalAlignment','center');
    posxlabel = get(xlabelh,'Position');
    posxlabel(1) = 0.5;
    set(xlabelh,'Position',posxlabel);
    
    set(ylabelh,'HorizontalAlignment','center');
    posylabel = get(ylabelh,'Position');
    posylabel(2) = 0.5;
    set(ylabelh,'Position',posylabel);
    
    set(gcf,'Units','centimeters');
    set(mainaxeshandles_,'Units','centimeters');
    set(xlabelh,'Units','centimeters');
    set(ylabelh,'Units','centimeters');
    
    % If the white space at the left and bottom are to be fixed, handle that
    if(fixedWhiteSpace_)
        
        oldPos_ = get(mainaxeshandles_,'Position');
        
        newPos_(3) = oldPos_(3) - boundaryWidth_(1) - boundaryWidth_(3);
        newPos_(4) = oldPos_(4) - boundaryWidth_(2) - boundaryWidth_(4);
        newPos_(1:2) = oldPos_(1:2) + boundaryWidth_(1:2);
        
        set(mainaxeshandles_,'Position',newPos_);
        
        set(gcf,'Units','normalized');
        set(mainaxeshandles_,'Units','normalized');
        set(xlabelh,'Units','normalized');
        set(ylabelh,'Units','normalized');
        
        set(xlabelh,'HorizontalAlignment','center');
        posxlabel = get(xlabelh,'Position');
        posxlabel(1) = 0.5;
        set(xlabelh,'Position',posxlabel);
        
        set(ylabelh,'HorizontalAlignment','center');
        posylabel = get(ylabelh,'Position');
        posylabel(2) = 0.5;
        set(ylabelh,'Position',posylabel);
        
        set(gcf,'Units','centimeters');
        set(mainaxeshandles_,'Units','centimeters');
        set(xlabelh,'Units','centimeters');
        set(ylabelh,'Units','centimeters');
        
    end
    
    set(gcf,'Units','normalized');
    set(xlabelh,'Units','normalized');
    set(ylabelh,'Units','normalized');
    set(mainaxeshandles_,'Units','normalized');
    set(legendaxeshandles_,'Units','normalized');
    
else
    %% Several axes, general adjustments
    % For several axes, meaning subplots, the task is more complicated.
    % In order to get good results, the user should specify the number of
    % subplots in x and y direction. Then they are aligned, and arranged so
    % they don't overlap.
    if(length(spreadSubplots_) ~= 2)
        % If nothing is specified, division by number of axes
        spreadSubplots_ = length(allaxeshandles_)*ones(1,2);
        %         spreadSubplots_ = length(allaxes_)*ones(1,2);
        if(verbose_ > 1)
            warning('Number of subplots not specified! You get better results using the option ''subplots'' together with the dimensions.');
        end
    end
    % Spread axes so they use all the space
    allOuterPos_ = zeros(length(mainaxeshandles_),4);
    for h_ = 1:length(mainaxeshandles_)
        oldUnits_ = get(mainaxeshandles_(h_),'Units');
        set(mainaxeshandles_(h_),'Units','normalized');
        OuterPos_ = get(mainaxeshandles_(h_),'OuterPosition');
        OuterPos_(1:2:3) = round(OuterPos_(1:2:3) * spreadSubplots_(2))/spreadSubplots_(2);
        OuterPos_(2:2:4) = round(OuterPos_(2:2:4) * spreadSubplots_(1))/spreadSubplots_(1);
        OuterPos_(3) = max(OuterPos_(3), 1/spreadSubplots_(2)); % Ensure positive width
        OuterPos_(4) = max(OuterPos_(4), 1/spreadSubplots_(1)); % Ensure positive height
        allOuterPos_(h_,:) = OuterPos_;
        set(mainaxeshandles_(h_),'OuterPosition',OuterPos_);
        set(mainaxeshandles_(h_),'Units',oldUnits_);
    end
    % Determine number of fields in x and y direction
    numX_ = 0;
    numY_ = 0;
    for g_ = 0:spreadSubplots_(2)
        numX_ = numX_ +any(any(allOuterPos_(:,1) == g_/spreadSubplots_(2)));
        numY_ = numY_ +any(any(allOuterPos_(:,2) == g_/spreadSubplots_(1)));
    end
    % Align subplots so all legends are visible
    allPosition_ = cell2mat(get(mainaxeshandles_,'Position'));
    for h_ = 0:numX_-1 % determine X-direction
        % Left boundary
        idx_ = allOuterPos_(:,1)==h_/spreadSubplots_(2);
        allPosition_(idx_,3) = allPosition_(idx_,3) - (max(allPosition_(idx_,1))-allPosition_(idx_,1));
        allPosition_(idx_,1) = max(allPosition_(idx_,1));
        % Right boundary
        idx_ = (allOuterPos_(:,1)+allOuterPos_(:,3))==(h_+1)/spreadSubplots_(2);
        minPos_ = min(allPosition_(idx_,1)+allPosition_(idx_,3));
        allPosition_(idx_,3) = minPos_-allPosition_(idx_,1);
    end
    for h_ = 0:numY_-1 % determine Y-direction
        % Lower boundary
        idx_ = allOuterPos_(:,2)==h_/spreadSubplots_(1);
        allPosition_(idx_,4) = allPosition_(idx_,4) - (max(allPosition_(idx_,2))-allPosition_(idx_,2));
        allPosition_(idx_,2) = max(allPosition_(idx_,2));
        % Upper boundary
        idx_ = (allOuterPos_(:,2)+allOuterPos_(:,4))==(h_+1)/spreadSubplots_(1);
        minPos_ = min(allPosition_(idx_,2)+allPosition_(idx_,4));
        allPosition_(idx_,4) = minPos_-allPosition_(idx_,2);
    end
    for g_ = 1:length(mainaxeshandles_) % Apply settings
        set(mainaxeshandles_(g_),'Position',allPosition_(g_,:));
    end
    % If the white space at the left and bottom are to be fixed, handle that
    if(fixedWhiteSpace_)
        if(numX_ == 1)
            for h_ = 1:length(mainaxeshandles_)
                % temporarily change units to those of the figure
                oldUnits_ = get(mainaxeshandles_(h_),'Units');
                set(mainaxeshandles_(h_),'Units',get(gcf,'Units'));
                newPos_ = get(mainaxeshandles_(h_),'Position');
                newPos_(3) = newPos_(3)+newPos_(1)-boundaryWidth_(1);
                newPos_(1) = boundaryWidth_(1);
                set(mainaxeshandles_(h_),'Position',newPos_);
                set(mainaxeshandles_(h_),'Units',oldUnits_);
            end
        end
        if(numY_ == 1)
            for h_ = 1:length(mainaxeshandles_)
                % temporarily change units to those of the figure
                oldUnits_ = get(mainaxeshandles_(h_),'Units');
                set(mainaxeshandles_(h_),'Units',get(gcf,'Units'));
                newPos_ = get(mainaxeshandles_(h_),'Position');
                newPos_(4) = newPos_(4)+newPos_(1:2)-boundaryWidth_(2);
                newPos_(2) = boundaryWidth_(2);
                set(mainaxeshandles_(h_),'Position',newPos_);
                set(mainaxeshandles_(h_),'Units',oldUnits_);
            end
        end
        % Check if this made problems
        oldUnits_ = get(mainaxeshandles_,'Units');
        set(mainaxeshandles_,'Units','normalized');
        if(any(any(cell2mat(get(mainaxeshandles_,'OuterPosition'))<-1e-4)))
            warning('Using the option ''fixedwhitespace'' causes bad results, try running %s without this option again!',mfilename);
        end
        for h_ = 1:length(mainaxeshandles_)
            set(mainaxeshandles_(h_),'Units',oldUnits_{h_});
        end
    end
end

% Set units of all axes to normalized, this enables the manual resizing
if newMatlab_
    set(mainaxeshandles_,'Units','normalized');
    set(legendaxeshandles_,'Units','normalized');
else
    set(allaxeshandles_,'Units','normalized');
end
% set(gcf,'Renderer',renderer_);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(length(figHandle_) > 1)
    % If several figures shall be treated, recursively call again
    idx_ = find(strcmpi('figure',varargin));
    myOpts = varargin;
    myOpts{idx_+1} = figHandle_(2:end);
    itm_formatfig(myOpts{:});
end

% set units for printing
set(gcf,'Units','Inches');
printpos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[printpos(3), printpos(4)])

% optional saving
if figname_ ~= 0
    if figpath_(end) ~= '/'
        figpath_ = [figpath_ '/'];
    end
    print(gcf,'-dpdf',[figpath_ figname_ '.pdf'],'-r0');
    %     print(gcf,'-depsc',[figpath_ figname_ '.eps'],'-r0');
    print(gcf,'-dpng',[figpath_ figname_ '.png'],'-r600');
    saveas(gcf,[figpath_ figname_ '.fig']);
    if(verbose_>0)
        fprintf('itm_formatfig: Figure saved\n');
    end
end


