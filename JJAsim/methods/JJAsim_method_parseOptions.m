function options = JJAsim_method_parseOptions(inputParameters,varInput,filename)
inputParameters = reshape(inputParameters',[],1);
options = struct(inputParameters{:});
optionNames = fieldnames(options);
nArgs = length(varInput);
if round(nArgs/2)~=nArgs/2
    error([filename,' propertyName/propertyValue pairs'])
end
for pair = reshape(varInput,2,[]) %# pair is {propName;propValue}
    inpName = pair{1};
    if any(strcmp(inpName,optionNames))
        options.(inpName) = pair{2};
    else
        error(['property ',inpName,' is not a valid option.'])
    end
end