function [raw_im,ff_im] = FlatFieldFilter(image_index)
if(image_index < 10^3)
    index_str = PadWithZeros(num2str(image_index),4);
else
    index_str = num2str(image_index);
end
folder = '/media/wez/Seagate Backup Plus Drive/nematic-liquid-crystals-8-30-17/images';
flat_field = double(imread(strcat(folder,'/Light.bmp')));
dark_frame = double(imread(strcat(folder,'/Dark.bmp')));
c = flat_field - dark_frame;
c = mean(mean(c))./c;
filename = strcat(folder,'/Image',index_str,'.bmp');
raw_im = double(imread(filename));
ff_im = c.*(raw_im-dark_frame);
end