function I = extractSpatialModes(im,window)
    [N,M] = size(im);
    nmin = window(1);
    nmax = window(2);
    mmin = window(3);
    mmax = window(4);
    W = zeros(N,M);
    W(nmin:nmax,mmin:mmax) = 1;
    F = fft2(im);
    WF = W.*F;
    I = ifft2(WF);
end