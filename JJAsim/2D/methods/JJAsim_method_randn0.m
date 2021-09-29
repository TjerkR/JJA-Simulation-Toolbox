function n = JJAsim_method_randn0(f,Np,W)
if size(f,1) > 1
    f = mean(f,1);
end
P = round(f*Np);
if length(P) == 1
    P = repmat(P,1,W);
end
n = zeros(Np,W,'int32');
for w = 1:W
    n(1:P(w),w) = 1;
    n(:,w) = n(randperm(Np),w);
end
end