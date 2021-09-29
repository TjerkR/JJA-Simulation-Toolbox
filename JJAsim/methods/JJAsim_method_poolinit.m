function JJAsim_method_poolinit(cores)
if cores > 1
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(cores)
    else
        if poolobj.NumWorkers ~= cores
            delete(gcp('nocreate'));
            parpool(cores);
        end
    end
end
end