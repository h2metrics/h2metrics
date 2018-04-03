%Input 'filename.txt' and 'figfile.bmp'
function ConvertTxtToBmp(loadfile,figfile)
    X=load(loadfile);
    fill(X(2:end,1),X(2:end,2),'k');
    saveas(gcf,figfile);
end