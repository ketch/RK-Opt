function wf=writeField(writeFid,name,value)
%function wf=writeField(name,value)

fprintf(writeFid,'\n%s\n',name);
[rows cols] = size(value);
x = repmat('%5.16E\t',1,(cols-1));
fprintf(writeFid,[x,'%5.16E\n'],value);

