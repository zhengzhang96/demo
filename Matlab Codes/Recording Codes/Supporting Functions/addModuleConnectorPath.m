function addModuleConnectorPath(optarg)
    addpath('../../Matlab Codes/');
    addpath('../../Libraries/include/');
    if ispc
        if nargin == 1
            if strcmp(optarg,'win32')
                addpath('../../Libraries/lib32/');
                return;
            end
        else
            addpath('../../Libraries/lib64/');
        end
    elseif ismac
        addpath('../../Libraries/libmac/');
    elseif islinux
        addpath('../../Libraries/liblinux/');
    end
end