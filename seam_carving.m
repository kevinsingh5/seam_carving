
% Load image and convert to double and grayscale
image=imread('castle.jpeg');
%image=im2double(image);
%image_gray=rgb2gray(image);

vert_image = removeVertical(image, 100);

% Compute the energy function at each pixel using the magnitude of the 
%   x and y gradients
%[fx, fy] = imgradientxy(image_gray, 'prewitt'); % equivalent to doing dF/dx and dF/dy
%imshowpair(fx, fy);
%energy_map = sqrt(fx.^2 + fy.^2); 
%imshow(energy_map);
%imagesc(energy_map);
%colorbar;

%new_image = removeVertical(image, 100);
%imshow(new_image);

