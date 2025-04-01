function h = my_figure(id)


fntSize = 14;
set(0,'defaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', fntSize);
set(0,'defaultTextFontName', 'Times New Roman')
set(0,'defaultTextFontSize', fntSize)


h0 = figure(id);
box on
set(gca,'fontsize',14)
set(gca,'fontname','Times New Roman')
% set(gca,'textfontsize',14)

if(nargout == 1)
    h = h0;
end


end