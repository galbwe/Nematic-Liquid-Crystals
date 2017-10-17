function cm = findCritModes(P,window)    
    N1 = window(1);
    N2 = window(2);
    M1 = window(3);
    M2 = window(4);
    [N,M] = size(P);
    W = zeros(N,M);
    W(N1:N2,M1:M2) = 1;
    PW = P.*W;
    NP = sum(sum(PW(N1:N2,M1:M2)));
    [m,n] = meshgrid(M1:M2,N1:N2);
    mc = round(sum(sum(m.*PW(N1:N2,M1:M2)))/NP);
    nc = round(sum(sum(n.*PW(N1:N2,M1:M2)))/NP);
    cm = [nc,mc];
end