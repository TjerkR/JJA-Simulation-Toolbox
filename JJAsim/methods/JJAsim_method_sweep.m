function xProblemList = JJAsim_method_sweep(Wx,xDim,Ws,rep)
sizes = [1,cumprod(Ws)];
step = sizes(xDim);
if Wx == 1
    probs = ones(1,Ws(xDim));
else
    if Ws(xDim) ~= Wx
       error('Wx must be 1 or Ws(xDim)');
    end
    probs = 1:Wx;
end
xProblemList = reshape(ones(step,1)*probs,[],1);
xProblemList = repmat(xProblemList,rep*prod(Ws(xDim+1:end)),1);
end