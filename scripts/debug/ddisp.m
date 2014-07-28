% Simple debug display function

function ddisp(varargin)
    
    global debug;
    
    if ~debug
        return
    end
    
    for i = 1:nargin
        disp(varargin{i});
    end
end
