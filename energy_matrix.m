%% This function computes the energy function on a given image and 
%   returns the resulting [double, grayscale] energy image

function energy_image = energy_matrix(image)
% Convert image to double and grayscale; add smooth filter
gray_image = rgb2gray(image);
image_double = im2double(image);
gray_double = rgb2gray(image_double);

% testing Guassian smoothing
hsize = 4;
sigma = 2;
h = fspecial('gaussian', hsize, sigma);
%mesh(h);
%imagesc(h);
smooth = imfilter(gray_image, h);
dbsmooth = im2double(smooth);
%imshow(dbsmooth);

% testing with convolution filter
fx = [-1, 1];
fy = [-1; 1];
gradx = double(imfilter(gray_image, fx, "conv", "replicate"));
grady = double(imfilter(gray_image, fy, "conv", "replicate"));
%imshowpair(gradx, grady);


%% Compute the energy function at each pixel using the magnitude of the x and y gradients
[gx, gy] = imgradientxy(gray_double);         % equivalent to doing dF/dx and dF/dy
%imshowpair(gx, gy);
%[energy_map, gdir] = imgradient(dbsmooth);   % another way to get energy map
energy_map = sqrt((gx.^2 + gy.^2)); 
%energy_map = abs(gx) + abs(gy);
%imshow(energy_map);

energy_image = energy_map;

end

