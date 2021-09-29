function xProblemList = JJAsim_method_list(Wx,W,rep)
if Wx == 1
    xProblemList = ones(W*rep,1);
else
    if Wx == W
        xProblemList = repmat((1:W)',rep,1);
    else
       error('Wx must be 1 or W'); 
    end
end