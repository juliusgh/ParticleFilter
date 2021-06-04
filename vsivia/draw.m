function draw(x, y, varargin)

    if numel(varargin) > 1

        arrayfun(@(xl,xu,yl,yu) rectangle('Position',[xl yl xu-xl yu-yl],'LineWidth',0.5,'LineStyle','-','EdgeColor',varargin{2},'FaceColor',varargin{1}), x.lower, x.upper, y.lower, y.upper) ;
        
    else
        
        %arrayfun(@(xl,xu,yl,yu) rectangle('Position',[xl yl xu-xl yu-yl],'LineWidth',0.5,'LineStyle','-','EdgeColor','k','FaceColor',varargin{1}), x.lower, x.upper, y.lower, y.upper) ;
        arrayfun(@(xl,xu,yl,yu) rectangle('Position',[xl yl xu-xl yu-yl],'LineWidth',0.5,'LineStyle','none','FaceColor',varargin{1}), x.lower, x.upper, y.lower, y.upper) ;

    end

end
