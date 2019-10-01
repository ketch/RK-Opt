function writeField(writeFid,name,value)
% function writeField(writeFid,name,value)

fprintf(writeFid,'\n%s\n',name);
fprintf(writeFid, [repmat('%5.16E\t', 1, size(value,2)),'\n'], value'); 

