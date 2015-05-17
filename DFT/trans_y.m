function new_Y = trans_y(Y)

% apply the DFT to each column of Y
new_Y = fft(Y)/sqrt(size(Y, 1));
    
end