%% plot hole locations over the pattern for several different threshold
%% values
clear
close all
clc
load('mean_power.mat')
mean_power = fftshift(mean_power2);
[N,M] = size(mean_power);
window = [454,467,14,36];%window for mode 2
cm = findCritModes(mean_power,window);
im_index = 2017;
[raw_im,ff_im] = FlatFieldFilter(im_index);
multiple_thresholds = figure(1);
set(multiple_thresholds,'position',[0,0,1600,800]);
thresholds = [.02,.04,.06,.08,.1,.12,.14,.16];
for i = 1:8
    thresh = thresholds(i);
    [I,A,H] = findEnvAndHoles(ff_im,window,cm,thresh);
    [X,Y] = findOnes(H);
    subplot(2,4,i);
    imagesc(real(I));
    colormap('hot');
    hold on
    scatter(X,Y,8,'g','filled')
    title(strcat('threshold = ',num2str(thresh)))
end
%% same thing, but as a video
clear
close all
clc
load('mean_power.mat')
mean_power = fftshift(mean_power2);
[N,M] = size(mean_power);
window = [454,467,14,36];%window for mode 2
cm = findCritModes(mean_power,window);
im_index = 10000;
[raw_im,ff_im] = FlatFieldFilter(im_index);
multiple_thresholds = figure(2);
set(multiple_thresholds,'position',[0,0,1600,800]);
thresholds = .001:.001:.1;
for i = 1:length(thresholds)
    thresh = thresholds(i);
    [I,A,H] = findEnvAndHoles(ff_im,window,cm,thresh);
    [X,Y] = findOnes(H);
    imagesc(real(I));
    colormap('hot');
    hold on
    scatter(X,Y,8,'g','filled')
    title(strcat('threshold = ',num2str(thresh)))
    hold off
    pause(1)
end
