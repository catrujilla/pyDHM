clc
clear all 
close all

%%
filename = 'I3_Fresnel_Lens.bmp';
holo = double(imread(filename));
holo = holo(:,:,1);
cutImage = uint8(255 * mat2gray(holo));   
imwrite(cutImage,gray(255),'gbCells.jpg','jpg')
[M,N] = size(holo);
if (M > N)
    cut = (M - N)/2;
    holo = holo(cut:(M-cut-1):N,1:N);
elseif (M < N)
    cut = (N - M)/2;
    holo = holo(1:M,cut:(N-cut)-1);
else
    holo = holo;
end

%% save
holo = uint8(255 * mat2gray(holo));
imwrite(holo,gray(255),filename,'jpg')