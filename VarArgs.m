classdef VarArgs
    properties
        argstruct = struct;
        ignorecase = 0;
    end
    
    methods
        function obj = VarArgs(vararginlist, varargin)
            % parse varargins of this function
            i = 1;
            while i < length(varargin)
                if ischar(varargin{i})
                    argname = lower(varargin{i});
                    switch argname
                        case 'ignorecase'
                            argval = TRUE;
                            if ~ischar(varargin{i+1})
                                i=i+1; argval = logical(varargin{i});
                            end
                            obj.ignorecase = argval;
                    end
                end
            end
            % parse the vararginlist given
            nargs = length(vararginlist);
            obj.argstruct = struct;
            argvals = cell(1,nargs);
            argname = 'defaults';
            nargvals = 0;
            for i = 1:nargs
                varval = vararginlist{i};
                if ischar(varval) && (nargvals > 0 || i == 1)
                    if nargvals == 1
                        obj.argstruct.(argname) = argvals{1};
                    elseif nargvals > 0
                        obj.argstruct.(argname) = argvals(1:nargvals);
                    end
                    nargvals = 0;
                    argname = varval;
                    if obj.ignorecase; argname = lower(argname); end
                else
                    nargvals = nargvals + 1;
                    argvals{nargvals} = varval;
                end
            end
            if nargvals == 1
                obj.argstruct.(argname) = argvals{1};
            elseif nargvals > 0
                obj.argstruct.(argname) = argvals(1:nargvals);
            end
        end
        
        function retval = get(obj, paramname, defaultval)
            if obj.ignorecase; paramname = lower(paramname); end
            if isfield(obj.argstruct, paramname)
                retval = obj.argstruct.(paramname);
            else
                retval = defaultval;
            end
        end
    end
    
end

