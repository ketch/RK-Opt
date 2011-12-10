function wfb=writeFile_butcher(A,b,c,class,writeFid)
%function wfb=writeFile_butcher(A,b,c,class,writeFid)
%
% 
% Write to file Butcher's coefficients

str = 'A';
if writeFid ~= -1
    fprintf(writeFid,'%s\r\n',str);
end
[rows cols] = size(A);
x = repmat('%5.16E\t',1,(cols-1));
fprintf(writeFid,[x,'%5.16E\n'],A');

str = 'b';
fprintf(writeFid,'\n%s\r\n',str);
[rows cols] = size(b');
x = repmat('%5.16E\t',1,(cols-1));
fprintf(writeFid,[x,'%5.16E\n'],b');

str = 'c^T';
fprintf(writeFid,'\n%s\r\n',str);
[rows cols] = size(c');
x = repmat('%5.16E\t',1,(cols-1));
fprintf(writeFid,[x,'%5.16E\n\n'],c');

if (class(1:2)=='2S' | class(1:2)=='3S')
    str = '==============================================================';
    fprintf(writeFid,'\n%s\r\n\n',str);
end

wfb = 1;

