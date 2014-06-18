% This is a simple MATLAB program to easily initialize a mutex lock
% that can be used to share data between threads using MATLAB
% parfor loops. This implementation is neither complex nor robust
% and is merely a convient abstraction that works MOST of the time.
%
% Written by Andrew Gilchrist-Scott

function mutexInit()
    global = lockDict;
    
    lockDict = containers.Map;
end

function mutexCreate(key)
    global lockDict;
    
    lockDict(key) = unit8(0)
end

function mutexGet(key)
    global lockDict;
    
    lockGotten = false;
    
    while ~lockGotten
        lockDict(key) = lockDict(key) + 1;
        % if another thread does this simultaneously, both are
        % locked; note that this will keep more than one thread
        % from gaining the mutex at the same time, but it will make
        % for a long wait if two do try to grab the lock at the
        % same time and no one has the lock
        if lockDict(key) == 1
            lockGotten = true;
        else
            pause(.02)
        end
    end
end

function mutexRelease(key)
    global lockDict;
    
    lockDict(key) = 0;
end

function mutexDestroy();
    global lockDict;
    
    clear lockDict
end


        