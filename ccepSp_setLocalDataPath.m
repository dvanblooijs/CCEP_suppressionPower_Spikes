function localDataPath = ccepSp_setLocalDataPath(varargin)
%
% input:  
%   ccepSp_personalDataPath: optional, set to 1 if adding ccepSp_personalDataPath
%
% when adding ccepSp_personalDataPath, the following function should be in the
% root of this repo:
%
% function localDataPath = ccepSp_personalDataPath()
%     'localDataPath = [/my/path/to/data];
%             
% ccepSp_personalDataPath is ignored in .gitignore
%

if ~isempty(varargin)
    % add path to data
    if varargin{1}==1 && exist('ccepSp_personalDataPath','file')
        localDataPath = ccepSp_personalDataPath();
                
    elseif varargin{1}==1 && ~exist('ccepSp_personalDataPath','file')
        sprintf(['add ccepSp_personalDataPath function to add your localDataPath:\n'...
            '\n'...
            'function localDataPath = ccepSp_personalDataPath()\n'...
            '\n'...
            'localDataPath.proj_dirinput = [/my/path/to/data];\n'...
            'localDataPath.proj_diroutput = [/my/path/to/output];\n'...
            '\n'...
            'this function is ignored in .gitignore'])
        return
    end
    
end

return

