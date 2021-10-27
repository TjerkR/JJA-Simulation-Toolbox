function Ic_array = simple_array2_to_Ic_array(array)
if size(array,1) ~= 2*size(array,2)
    error("Array is not the right size.")
end
N = size(array,2);

% 1
working_array1 = flipud(array);

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