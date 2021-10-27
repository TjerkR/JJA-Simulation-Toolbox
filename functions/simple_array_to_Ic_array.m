function Ic_array = simple_array_to_Ic_array(array)
if size(array,1) ~= size(array,2)
    error("Array is not square.")
end
N = length(array);

% 0
working_array0 = flipud(array);

% 1
working_array1 = zeros(N*2, N);
for i = 1:(N*2)
    for j = 1:N
        working_array1(i,j) = working_array0(ceil(i/2),j);
    end
end

% 2
working_array2 = zeros(N*2, N+1);
for i = 1:(N*2)
    if mod(i,2) ~= 0
        working_array2(i,:) = [0 working_array1(i,:)];
    else
        working_array2(i,:) = [working_array1(i,:) 0];
    end
end

% 3
working_array3 = working_array2(1:(N*2-1), 1:N);

% final
Ic_array = flipud(working_array3);

end