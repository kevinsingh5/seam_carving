% This function removes vertical seams from an image

function output_image = removeVertical(im, numPixels)
%% Convert image to double and grayscale; add smooth filter
gray = rgb2gray(im);
image = im2double(im);
gray_image = rgb2gray(image);

% testing Guassian smoothing
hsize = 4;
sigma = 2;
h = fspecial('gaussian', hsize, sigma);
%mesh(h);
%imagesc(h);
smooth = imfilter(gray, h);
dbsmooth = im2double(smooth);
%imshow(dbsmooth);

% testing with convolution filter
fx = [-1, 1];
fy = [-1; 1];
gradx = double(imfilter(gray, fx, "conv", "replicate"));
grady = double(imfilter(gray, fy, "conv", "replicate"));
%imshowpair(gradx, grady);


%% Compute the energy function at each pixel using the magnitude of the x and y gradients
[gx, gy] = imgradientxy(gray_image);    % equivalent to doing dF/dx and dF/dy
%imshowpair(gx, gy);
%[energy_map, gdir] = imgradient(dbsmooth);   % another way to get energy map
energy_map = sqrt((gx.^2 + gy.^2)); 
%energy_map = abs(gx) + abs(gy);
%imshow(energy_map);


%% Create a new matrix to store min energy seams
[rows, cols] = size(energy_map);
M = zeros(rows, cols);  % new matrix for cumuluative min energy seams
M(:) = energy_map(:);

for j=1:cols   % iterate cols
    for i=2:rows   % iterate rows
        if j==1
            M(i, j) = energy_map(i, j) + min([M(i-1, j), M(i-1, j+1)]);
            continue;
        elseif j==cols
            M(i, j) = energy_map(i, j) + min([M(i-1, j-1), M(i-1, j)]);
            continue;
        end
        M(i, j) = energy_map(i, j) + min([M(i-1, j-1), M(i-1, j), M(i-1, j+1)]);
    end
end

% for i=2:rows
%    for j=2:cols-1
%        %if(j==1)
%        %    M(i, j) = energy_map(i, j) + min([M(i-1, j), M(i-1, j+1)]);
%        %    continue;
%        %elseif(j==cols)
%        %    M(i, j) = energy_map(i, j) + min([M(i-1, j-1), M(i-1, j)]);
%        %    continue;
%        %else
%        M(i, j) = energy_map(i, j) + min([M(i-1, j-1), M(i-1, j), M(i-1, j+1)]);
%        %end
%    end
%    M(i,1) = energy_map(i,1) + min([M(i-1,1), M(i-1,2)]);
%    M(i,cols) = energy_map(i,cols) + min([M(i-1,cols), M(i-1,cols-1)]);
% end

%imshow(M);


%% get min seam in energy map 
new_image = zeros(size(image));
new_image(:) = image(:);

for i=1:numPixels   % decrease the image by the given # of pixels
%     minPx = M(rows, 1);
%     minCol = 1;
%     for j=1:cols     % find the min value in last row of M
%         if(M(rows, j) < minPx)
%             %minPx = (M(rows, j));
%             %minCol = j;     
%             M(rows, minCol) = 1.00;    % take that px out of consideration
%         end
%     end
    [minPx, minCol] = min(M(rows,:));
    M(rows, minCol) = 1.0;    % take that px out of consideration

    % trace back the min seam to remove
    for k=rows-1:-1:1
        if((minCol == 1) && M(k, minCol) <= M(k, minCol+1))
            minCol = minCol;
            continue;
        elseif((minCol == 1) && M(k, minCol+1) <= M(k, minCol))
            minCol = minCol+1;
            continue;
        elseif(minCol == cols && M(k, minCol) <= M(k, minCol-1))
            minCol = minCol;
            continue;
        elseif(minCol == cols && M(k, minCol-1) <= M(k, minCol))
            minCol = minCol-1;
            continue;
        elseif(M(k, minCol+1) < M(k, minCol-1) && M(k, minCol+1) < M(k, minCol))  % check if top right px is the minimum one
            minCol = minCol+1;
            continue;
        elseif(M(k, minCol-1) < M(k, minCol) && M(k, minCol-1) < M(k, minCol+1)) % check if top left px is min one
            minCol = minCol-1;
        else % top px is the min one
            minCol = minCol;
        end
        M(k+1, minCol) = 1.0; % remove the min px from top row
        new_image(k+1, minCol) = 1.0;   % remove comment to show seams on image
    end
end

imshow(new_image);

output_image = new_image;

end

