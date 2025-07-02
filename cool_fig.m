function cool_fig(h)
    if isempty(findobj('Number',h))
        figure(h);
    else
        set(groot,'CurrentFigure',h);
    end
end