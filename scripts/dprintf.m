% this is a debug print statement to be used in all our code

function dprintf(varargin)
    
    global debug;
    
    if ~debug
        return
    end
    
    printArgs = zeros(1,nargin-1);
    
    for i = 1:nargin-1
        printArgs(i) = varargin{i+1};
    end
    
    fprintf(varargin{1},printArgs);
end

    