%% This function computes the cummulative minimun energy map for the
%       given energy map

function map = cmin_energy(energy_image)
%% Create a new matrix (energy map) to store min energy seams
[rows, cols] = size(energy_image);
M = zeros(rows, cols);  % new matrix for cumuluative energy map
M(:) = energy_image(:);

for j=1:cols   % iterate cols
    for i=2:rows   % iterate rows
        if j==1
            M(i, j) = energy_image(i, j) + min([M(i-1, j), M(i-1, j+1)]);
            continue;
        elseif j==cols
            M(i, j) = energy_image(i, j) + min([M(i-1, j-1), M(i-1, j)]);
            continue;
        end
        M(i, j) = energy_image(i, j) + min([M(i-1, j-1), M(i-1, j), M(i-1, j+1)]);
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

map = M;

end

