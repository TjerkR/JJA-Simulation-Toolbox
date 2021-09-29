function xP = JJAsim_method_tocell(x,xTableSelection,wDim,nrOfDimensions)
if wDim > nrOfDimensions
    error('working dimension exceeds the number of dimensions');
end
if nrOfDimensions > 6
    error('number of dimensions larger than 6, not supported');
end
partitions = size(xTableSelection,2);
xP = cell(1,partitions);
for i = 1:partitions
    switch wDim
        case 1
            xP{i} = x(xTableSelection(:,i),:,:,:,:,:);
        case 2
            xP{i} = x(:,xTableSelection(:,i),:,:,:,:);
        case 3
            xP{i} = x(:,:,xTableSelection(:,i),:,:,:);
        case 4
            xP{i} = x(:,:,:,xTableSelection(:,i),:,:);
        case 5
            xP{i} = x(:,:,:,:,xTableSelection(:,i),:);
        case 6
            xP{i} = x(:,:,:,:,:,xTableSelection(:,i));
    end
end
end