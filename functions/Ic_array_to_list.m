function Ic = Ic_array_to_list(array)
N = size(array, 2);
Nj = 2*N*(N-1);

rows = N*2-1;
cols = N;

Ic = ones(Nj, 1);
index_array = Ic_list_to_array(1:Nj, 0);

for i = 1:rows
    for j = 1:cols
        if index_array(i,j) ~= 0
            Ic(index_array(i,j)) = array(i,j);
        end
    end
end

end