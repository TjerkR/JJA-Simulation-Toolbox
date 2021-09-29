function [x,xCompact] = JJAsim_method_checkInput(x,xType,xSize,xSizeType,xName)
%a dimension is compact if its size is 1 despite its xSize being larger than 1. 

%default values for optional input
if nargin < 4
    xSizeType = zeros(size(xSize));
end
if nargin < 5
    xName = 'variable';
end

%check for nan or inf values
if sum(reshape(isnan(x),[],1)) ~=0
    error([xName,' contains NaN entries']);
end
if sum(reshape(isinf(x),[],1)) ~=0
    error([xName,' contains inf entries']);
end

%cast to required class if necessary
xClass = class(x);
if strcmp(xType,'matlabInt')
    if ~strcmp(xClass,'double') %#ok<STISA>
        if ischar(x)
            error([xName,' is a char whereas ',xType,' is required.']);
        end
        warning(['converting ',xName,' from ',xClass,' to double.']);
        x = cast(x,xType);
    end
    if ~isequal(x,round(x)) || sum(reshape(x <= 0,[],1))~=0
       error(['variable ',xName,' must contain positive whole numbers']) 
    end
else
    if ~strcmp(xClass,xType)
        if ischar(x) && ~strcmp(xType,'char')
            error([xName,' is a char whereas ',xType,' is required.']);
        end
        warning(['converting ',xName,' from ',xClass,' to ',xType,'.']);
        x = cast(x,xType);
    end
end
xDims = length(xSize);

%get number of dimensions
nDims = ndims(x);
if nDims == 2 && size(x,2) == 1
    nDims = 1;
end

%check if x does not have too many dimensions
if nDims > xDims
    error([xName,' has too many dimensions.']);
end

%check the size of each dimension and if dimensions are compact. 
xCompact = false(xDims,1);
compactify = false(xDims,1);
xscale = mean(mean(mean(mean(mean(mean(mean(abs(x))))))));
if xscale == 0
    xscale = 1;
end
if isa(x,'single')
    xscale = xscale*1E-6;
else
    xscale = xscale*1E-12;
end

for xDim = 1:xDims
   switch xSizeType(xDim)
       case 0
            if size(x,xDim) ~= xSize(xDim)
                error([xName,' size in dimension ',num2str(xDim),' does not equal the required size.'])
            end
       case 1
            if ~(size(x,xDim) == xSize(xDim) || size(x,xDim) == 1)
                error([xName,' size in dimension ',num2str(xDim),' does not equal the required size.'])
            end
       case 2
           
       otherwise
           error('unrecognized xSizeType');
   end
   
   %set if x along xDim is compact
   if size(x,xDim) == 1
       xCompact(xDim) = true;
   end
   
   %check if x can be compactified along dimension
   if size(x,xDim) > 1 && xSizeType(xDim) ~= 0
       pdif = mean(mean(mean(mean(mean(mean(abs(std(x,[],xDim))))))));
      if pdif < xscale
          compactify(xDim) = true;
          xCompact(xDim) = true;
      end
   end
end

%compactify x
scell = cell(1,xDims);
for xDim = 1:xDims
    if compactify(xDim)
        scell{xDim} = 1;
    else
        scell{xDim} = 1:size(x,xDim);
    end
end
x = x(scell{:});
end