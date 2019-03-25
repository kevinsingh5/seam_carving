function new_image = removeHorizontal(image,cut_cols)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[rows, cols] = size(image);
new_cols = cols - cut_cols;
M = zeros(rows, cols);

M(:,1) = image(:,1);
for j=2 : cols-1
    for i=2 : rows-1   
        M(i,j) = image(i,j) + min([M(i-1,j-1), M(i,j-1), M(i+1,j-1)]);
    end
    M(1,j) = image(1,j) + min([M(1,j-1), M(2,j-1)]);
    M(rows,j) = image(rows,j) + min([M(rows-1,j-1), M(rows,j-1)]);
end

new_image = M;
%imshow(new_image);
end

