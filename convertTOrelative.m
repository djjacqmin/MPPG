function rel = convertTOrelative(x,y,xnorm)

if xnorm == 'max'
    ynorm = max(y);
else
    ynorm = interp1(x,y,xnorm);
end

rel = y ./ ynorm * 100;