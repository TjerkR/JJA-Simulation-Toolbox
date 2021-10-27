function [Ic_list, array, sigma_in_actual, mu_in, mu_in_actual] = generate_disordered_Ic(N, sigma, NameValueArgs)
arguments
    N (1,1) double
    sigma (1,1) double
    NameValueArgs.mu = 1
    NameValueArgs.l_disorder = 0
    NameValueArgs.sigma_in (1,1) double = 0
end

mu = NameValueArgs.mu;
l_disorder = NameValueArgs.l_disorder;
sigma_in = NameValueArgs.sigma_in;

if l_disorder == 0
    if sigma_in ~= 0
        warning("sigma_in is not used when l_disorder == 0.")
    end
    
    working_array = sigma * randn(N*2,N) + mu;
    
    sigma_in_actual = sqrt(var(working_array));
    mu_in_actual = mean(working_array);
    mu_in = mu;
else
    working_array = zeros(N*2,N);
    n = 0;
    for i = 1:2:(N*2-1)
        for j = 1:N
            if mod(ceil(i/2)+l_disorder-1, l_disorder) == 0 ...
                    && mod(j+l_disorder-1, l_disorder) == 0
                % Generate average Ic for this area
                clear Ic
                m = 0;
                n = n + 1;
                if l_disorder < N
                    mu_in(n) = sigma * randn() + mu;
                else
                    mu_in(n) = mu;
                end
                for ni = 1:(l_disorder*2)
                    for nj = 1:l_disorder
                        if ~(i+ni-1 > N*2 || j+nj-1 > N)
                            % Generate Ic for this junction
                            m = m + 1;
                            Ic(m) = sigma_in * randn() + mu_in(n);
                            working_array(i+ni-1, j+nj-1) = Ic(m);
                        end
                    end
                end
                sigma_in_actual(n) = sqrt(var(Ic));
                mu_in_actual(n) = mean(Ic);
            end
        end
    end
end

array = flipud(working_array);
Ic_array = simple_array2_to_Ic_array(array);
Ic_list = Ic_array_to_list(Ic_array);

end

