function padded_str = PadWithZeros(str,len)
    padded_str = str;
    while length(padded_str) < len
       padded_str = strcat('0',padded_str);
    end
end