%given values of spatial fourier transform inside a window and critical
%modes at center of fourier spectrum, creates an N by M array that is zero
%outside of the given window and takes the given values inside the window
function FS = fourierSpectrum(window,windowValues,spectrumSize,critModes)
    FS = zeros(spectrumSize);
    FS(window(1):window(2),window(3):window(4)) = windowValues;
    N = spectrumSize(1); M = spectrumSize(2);
    wn1 = window(1); wn2 = window(2);
    wm1 = window(3); wn2 = window(4);
    cn = critModes(1); cm = critModes(2);
    FS = FS([cn:N 1:cn-1],[cm:M 1:cm-1]);    
end