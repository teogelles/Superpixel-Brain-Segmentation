function [NLL,g] = UGM_CRF_NLL(w,Xnode,Xedge,Y,nodeMap,edgeMap,edgeStruct,inferFunc,varargin)

[nNodes,maxState] = size(nodeMap);
nNodeFeatures = size(Xnode,2);
nEdgeFeatures = size(Xedge,2);
nEdges = edgeStruct.nEdges;
edgeEnds = edgeStruct.edgeEnds;
nStates = edgeStruct.nStates;

nInstances = size(Y,1);
NLL = 0;
g = zeros(size(w));

for i = 1:nInstances
    
    % Make potentials
    if edgeStruct.useMex
        [nodePot,edgePot] = UGM_CRF_makePotentialsC(w,Xnode,Xedge,nodeMap,edgeMap,nStates,edgeEnds,int32(i));
    else
        [nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,i);
    end
    
    %     ignoreStd = std(reshape(edgePot,1,numel(edgePot))) + 1;
    %     ignoreMean = mean(mean(mean(edgePot)));
    %     edgePot = (edgePot-ignoreMean)/ignoreStd ;
    %     edgePot = edgePot + abs(min(min(min(edgePot)))) + 1;
    %     
    %     ignoreStd = std(reshape(nodePot, 1, numel(nodePot))) + 1;
    %     ignoreMean = mean(mean(mean(nodePot)));
    %     nodePot = (nodePot-ignoreMean)/ignoreStd ;
    %     nodePot = nodePot + abs(min(min(min(nodePot)))) + 1;
    
    % Compute marginals and logZ
    % if (isnan(nodePot(1, 1)))
    %     disp('nodePot Bad')
    %     disp(nodePot)
    %     disp(edgeStruct)
    %     disp(varargin)
    % end
    % if (isnan(edgePot(1, 1, 1)))
    %     disp('edgePot Bad')
    %     disp(edgePot)
    %     disp(edgeStruct)
    %     disp(varargin)
    % end

    disp('variable status pre-inferFunc')
    if (isnan(nodePot(1, 1)))
        disp('nodePot isnan')
    end
    if (isnan(edgePot(1, 1)))
        disp('edgePot isnan')
    end
    if (isnan(edgeStruct.V(1)))
        disp('edgeStruct V isnan')
    end
    if (isnan(edgeStruct.E(1)))
        disp('edgeStruct E isnan')
    end
    if (isnan(edgeStruct.edgeEnds(1, 1)))
        disp('edgeStruct edgeEnds isnan')
    end
    if (isnan(edgeStruct.nStates(1)))
        disp('edgeStruct nStates isnan')
    end
    [nodeBel,edgeBel,logZ] = inferFunc(nodePot,edgePot,edgeStruct,varargin{:});
    disp('variable status post-inferFunc')
    if (isnan(nodePot(1, 1)))
        disp('nodePot isnan')
    end
    if (isnan(edgePot(1, 1)))
        disp('edgePot isnan')
    end
    if (isnan(edgeStruct.V(1)))
        disp('edgeStruct V isnan')
    end
    if (isnan(edgeStruct.E(1)))
        disp('edgeStruct E isnan')
    end
    if (isnan(edgeStruct.edgeEnds(1, 1)))
        disp('edgeStruct edgeEnds isnan')
    end
    if (isnan(edgeStruct.nStates(1)))
        disp('edgeStruct nStates isnan')
    end
    
    % Update NLL
    if edgeStruct.useMex
        
        fprintf('\nNLL: %d, logZ: %d\n', NLL, logZ);
        if (isnan(logZ))
            disp('BadBadBad Logz')
        end
        
        NLL = NLL - UGM_LogConfigurationPotentialC(Y(i,:),nodePot, ...
                                                   edgePot,edgeEnds) + logZ;
        fprintf('\nNLL: %d, logZ: %d\n', NLL, logZ);
        
        % Updates in-place
        UGM_CRF_NLLC(g,int32(i),nodeBel,edgeBel,edgeEnds,nStates,nodeMap,edgeMap,Xnode,Xedge,Y);
        
        ignoreStd = std(reshape(g,1,numel(g))) + 1; %ADDED
        ignoreMean = mean(mean(mean(g)));
        g = (g-ignoreMean)/ignoreStd ;
        
    else
        NLL = NLL - UGM_LogConfigurationPotential(Y(i,:),nodePot,edgePot,edgeEnds) + logZ;
        
        if nargout > 1
            for n = 1:nNodes
                for s = 1:nStates(n)
                    for f = 1:nNodeFeatures
                        if nodeMap(n,s,f) > 0
                            if s == Y(i,n)
                                obs = 1;
                            else
                                obs = 0;
                            end
                            g(nodeMap(n,s,f)) = g(nodeMap(n,s,f)) + Xnode(i,f,n)*(nodeBel(n,s) - obs);
                        end
                    end
                end
            end
            for e = 1:nEdges
                n1 = edgeEnds(e,1);
                n2 = edgeEnds(e,2);
                for s1 = 1:nStates(n1)
                    for s2 = 1:nStates(n2)
                        for f = 1:nEdgeFeatures
                            if edgeMap(s1,s2,e,f) > 0
                                if s1 == Y(i,n1) && s2 == Y(i,n2)
                                    obs = 1;
                                else
                                    obs = 0;
                                end
                                g(edgeMap(s1,s2,e,f)) = g(edgeMap(s1,s2,e,f)) + Xedge(i,f,e)*(edgeBel(s1,s2,e) - obs);
                            end
                        end
                    end
                end
            end
        end
    end
end