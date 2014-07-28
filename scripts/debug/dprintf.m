% this is a debug print statement to be used in all our code

function dprintf(varargin)
    
    global debug;
    
    if ~debug
        return
    end
    
    printArgs = zeros(1,nargin-1);
    
    for i = 1:nargin-1
        if ~isstr(varargin{i+1})
            printArgs(i) = varargin{i+1};
        else
            %breaks if string isn't last item in list
            printArgs(i:i+length(varargin{i+1})-1) = varargin{i+1};
        end
    end
    
    fprintf(varargin{1},printArgs);
end

    