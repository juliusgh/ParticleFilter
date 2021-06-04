function plotIntervals(S,C,fig)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if nargin < 3
        fig = NaN;
    end
    if size(S,2) == 4
        plotIntervals(S(:,1:2),C,fig);
        plotIntervals(S(:,3:4),C,fig+1);
    else
        SB = [S.lower S.upper];
        SBu = unique(SB,'rows');
        X = [SBu(:,1)'; SBu(:,3)'; SBu(:,3)'; SBu(:,1)'];
        Y = [SBu(:,2)'; SBu(:,2)'; SBu(:,4)'; SBu(:,4)'];
        if isnan(fig)
            figure;
        elseif fig > 0
            figure(fig);
        end
        C_edge = C; %'k';
        patch(X,Y,C,'EdgeColor',C_edge,'FaceAlpha',1);
    end
end

