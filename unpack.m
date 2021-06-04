function varargout = unpack(x)
%UNPACK helper function to unpack intervals or arrays in the same way
%   The vsivia implementation for intervals uses a different indexing order
if isa(x,'interval')
    for k = 1:nargout
    	varargout{k} = x(:,k);
    end
else
    for k = 1:nargout
    	varargout{k} = x(k,:);
    end
end

