%% Apply flat field filter to raw data
M = 667;
N = 903;
T = 30000;
flat_field = double(imread('../raw-data/Light.bmp'));
dark_frame = double(imread('../raw-data/Dark.bmp'));
c = flat_field - dark_frame;
c = mean(mean(c))./c;
for i = 1:T
   index_str = num2str(i);
   while length(index_str) < 4
      index_str = strcat('0',index_str); 
   end
   raw_image_filename = strcat('../raw-data/Image',index_str,'.bmp');
   raw_image = double(imread(raw_image_filename));
   flat_field_image = c.*(raw_image-dark_frame);
   while length(index_str) < 5
      index_str = strcat('0',index_str); 
   end
   flat_field_image_filename = strcat('../flat-fielded-data/Image' ...
        ,index_str,'.bmp');
   imagesc(flat_field_image);
   colormap(gray)
   axis off
   saveas(gcf,flat_field_image_filename)
end
%% Read flat-fielded images into a .mat file
M = 667
N = 903
T = 30000
t_init = 1;
t_final = T;
flat_field_data = zeros(M,N,T);
for i = t_init:t_final
   i_str = PadWithZeros(num2str(i),5);
   filename = strcat('../flat-fielded-data/Image' ...
        ,i_str,'.bmp'); 
end

%% Take Fourier Transforms of flat-fielded images
M = 480;
N = 640;
T = 30000;
t_init = 1;
t_final = T;
fig1 = figure
fig2 = figure
fig3 = figure
for i = t_init:t_final
    i_str = PadWithZeros(num2str(i),5);
    filename = strcat('../flat-fielded-data/Image' ...
        ,i_str,'.bmp');
    flat_field = double(rgb2gray(imread(filename)));
    I = flat_field - mean(mean(flat_field));
    %figure(fig1)
    %imagesc(I);
    %colormap(gray)
    %axis off
    F = fft2(I);
    P = abs(fftshift(F))*2;%power spectrum
    figure(fig2)
    imagesc(P);
    axis off
    filename = strcat('../fourier-transform-data/Image' ...
        ,i_str,'.bmp');
    saveas(fig2,filename)    
end
%% compute 