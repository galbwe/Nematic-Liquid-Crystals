%assumed window = [nmin,nmax,mmin,mmax]
%assumed critModes = [nc,mc]
function [I,A,H] = findEnvAndHoles(im,window,critModes,thresh)
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
    nc = critModes(1);
    mc = critModes(2);
    A = ifft2(WF([nc:N,1:nc-1],[mc:M,1:mc-1]));
    H = (abs(A) < thresh);
    
end