function psd = pgm(x)
    N = length(x);
    psd = (1/N)*(abs(fft(x))).^2;
end