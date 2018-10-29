function [ varargout ] = dealarray( inarray )
n = length(inarray);
varargout = cell(1,n);
for i = 1:n
    varargout{i} = inarray(i);
end
end

