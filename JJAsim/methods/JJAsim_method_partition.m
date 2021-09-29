function [xProblemSelection,xTableSelection,Wx,partitionSize] = JJAsim_method_partition(...
    xProblemList,computePartitions)
%efficiently distributes problems over several partitions.
%example: xProblemList = [1 1 2 3 2 1 1]'
%         computePartitions = 4
%gives:   xProblemSelection = [1 1 2 1
%                              1 2 1 1]
%         xTableSelection = [1 0 1 1
%                            0 1 1 0
%                            0 1 0 0]
%         Wx = [1 2 2 1]
%         partitionSize = 2
%
% goes from (xProblemList,x,max(xProblemList)) to
% for p = 1:computePartitions
%       (xProblemSelection(:,p),x(xTableSelection(:,p),:),Wx(p))
% end
W = length(xProblemList);
WxMax = max(xProblemList);
partitionSize = ceil(W/computePartitions);
Wfull = partitionSize*computePartitions;
xProblemList = [xProblemList;repmat(xProblemList(end),Wfull-W,1)];
xProblemList = reshape(xProblemList,partitionSize,computePartitions);
xTableSelection = false(WxMax,computePartitions);
for w = 1:partitionSize
    for p = 1:computePartitions
        xTableSelection(xProblemList(w,p),p) = true;
    end
end
Wx = sum(xTableSelection,1);
xProblemSelection = zeros(size(xProblemList));
for  p = 1:computePartitions
    idx = find(xTableSelection(:,p));
    for w = 1:partitionSize
        xProblemSelection(w,p) = find(idx == xProblemList(w,p));
    end
end
end

