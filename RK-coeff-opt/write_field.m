function write_field(writeFid,name,value)
% function write_field(writeFid,name,value)
%
%
% Utility function to write a single parameter and value.

fprintf(writeFid,'\n%s\n',name);
fprintf(writeFid, [repmat('%5.16E\t', 1, size(value,2)),'\n'], value'); 
