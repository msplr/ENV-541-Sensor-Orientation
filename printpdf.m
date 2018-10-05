function printpdf(h,outfilename,xscale,yscale)
if nargin < 4
    yscale = 1;
end
if nargin < 3
    xscale = 1;
end
set(h, 'PaperUnits','centimeters');
set(h, 'Units','centimeters');
pos=get(h,'Position');
set(h, 'PaperSize', [xscale*pos(3) yscale*pos(4)]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition',[0 0 xscale*pos(3) yscale*pos(4)]);
print(outfilename, '-dpdf');
