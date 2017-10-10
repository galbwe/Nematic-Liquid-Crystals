function [peaks,i,j] = findpeaks2d(uhat,num_peaks,thresh,step)
    n0 = -2;
    n1 = -1;
    n2 = 0;
    count = 0;
    while n2 ~= num_peaks & count < 1000
        mask = (uhat > thresh);
        n2 = sum(sum(mask));
        if abs(n2 - n1) > 0 & n2 == n0
           step = step/2; 
        end
        if n2 > num_peaks
            thresh = thresh + step;
        elseif n2 < num_peaks 
            thresh = thresh - step;    
        end
        count = count + 1;
        n0 = n1;
        n1 = n2;        
    end
    peaks = uhat(uhat > thresh);
    i = zeros(1,length(peaks));
    j = zeros(1,length(peaks));
    for k = 1:length(peaks)
       [i0,j0] = find(uhat == peaks(k))
       i(k) = i0;
       j(k) = j0;
    end
end