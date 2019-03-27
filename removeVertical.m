% This function removes vertical seams from an image

function output_image = removeVertical(im, numPixels)
%% Convert image to double and grayscale
image = im2double(im);
gray_image = rgb2gray(image);

%% Compute the energy function at each pixel using the magnitude of the x and y gradients
%[gx, gy] = imgradientxy(gray_image); % equivalent to doing dF/dx and dF/dy
[energy_map, gdir] = imgradient(gray_image);
%imshowpair(gx, gy);
%imshow(energy_map);
%energy_map = sqrt((gx.^2 + gy.^2)); 
%energy_map = abs(gx) + abs(gy);
%imshow(energy_map);
%imagesc(energy_map);

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
%    for j=1:cols
%        if j==1
%            M(i, j) = energy_map(i, j) + min([M(i-1, j), M(i-1, j+1)]);
%            continue;
%        elseif j==cols
%            M(i, j) = energy_map(i, j) + min([M(i-1, j-1), M(i-1, j)]);
%            continue;
%        else
%            M(i, j) = energy_map(i, j) + min([M(i-1, j-1), M(i-1, j), M(i-1, j+1)]);
%        end
%    end
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
%             [minPx, minCol] = min(M(rows,:));
%             M(rows, minCol) = 1.00;    % take that px out of consideration
%         end
%     end
    [minPx, minCol] = min(M(rows,:));
    M(rows, minCol) = 1.0;    % take that px out of consideration
    % trace back the min seam to remove
    for k=rows-1:-1:1
        if(minCol == 1 && M(k, minCol) < M(k, minCol+1))
            minCol = minCol;
            continue;
        elseif(minCol == 1 && M(k, minCol+1) < M(k, minCol))
            minCol = minCol+1;
            continue;
        elseif(minCol == cols && M(k, minCol) < M(k, minCol-1))
            minCol = minCol;
            continue;
        elseif(minCol == cols && M(k, minCol-1) < M(k, minCol))
            minCol = minCol-1;
            continue;
        elseif(M(k, minCol+1) < M(k, minCol-1) && M(k, minCol+1) < M(k, minCol))  % check if top right px is the minimum one
            minCol = minCol+1;
        elseif(M(k, minCol-1) < M(k, minCol) && M(k, minCol-1) < M(k, minCol+1)) % check if top left px is min one
            minCol = minCol-1;
        else % top px is the min one
            minCol = minCol;
        end
        new_image(k+1, minCol) = 1.0; % remove the min px from top row
        M(k+1, minCol) = 1.0;
    end
end

imshow(new_image);

output_image = M;

end

