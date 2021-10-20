function array = Ic_list_to_array(list, L)
Nj = length(list);
N = 0.5 * (sqrt(2*Nj+1) + 1);

rows = N*2-1;
cols = N;
index_array = ones(rows, cols);
working_array = ones(rows, cols);

% creating index array
for i = 1:rows
    for j = 1:cols
        if mod(i,2) ~= 0 && j == 1
            index_array(i,j) = 0;
        elseif i == 1
            index_array(i,j) = 2*(j-1) - 1;
        elseif i == 2 && j ~= cols
            index_array(i,j) = 2*j;
        elseif i == 2 && j == cols
            index_array(i,j) = index_array(i,j-1) + 1;
        elseif mod(i,2) == 0
            index_array(i,j) = index_array(i-2, j) + rows;
        elseif mod(i,2) ~= 0 && i ~= rows
            index_array(i,j) = index_array(i-2, j) + rows;
        elseif i == rows && j == 2
            index_array(i,j) = index_array(i-2, j) + rows;
        elseif i == rows
            index_array(i,j) = index_array(i, j-1) + 1;
        end
    end
end

% creating Ic array from indices
for i = 1:rows
    for j = 1:cols
        if index_array(i,j) == 0
            working_array(i,j) = 0;
        else
            working_array(i,j) = list(index_array(i,j));
        end
    end
end

% removing hole
if L > 0
    hole_edge = (N-L)/2+1;
    working_array((hole_edge+0):(rows+1-hole_edge), hole_edge:(cols+1-hole_edge)) = 0;
    for i = 1:rows
        if mod(i,2) ~= 0 && i < (rows+1-hole_edge) && i > hole_edge
            working_array(i, cols+2-hole_edge) = 0;
        end
    end
end

% flipping to get proper array
array = flipud(working_array);

end