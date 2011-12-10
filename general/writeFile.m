function wf=writeFile(rk,p,class)
%function wf=writeFile(rk,p,ls)
%
% 
% Write to file Butcher's coefficients and low-storage coefficients if 
% required.

szA = size(rk.A);

outputFileName = strcat('ERK-',num2str(p),'-',num2str(szA(1)),'.txt');
writeFid = fopen(outputFileName,'w');

fprintf(writeFid, '%s\t\t %s\n', '#stage','order');
output = [szA(1);p];
fprintf(writeFid, '%u\t \t\t%u\n\n',output);

% Butcher coefficients
str = 'A';
if writeFid ~= -1
    fprintf(writeFid,'%s\r\n',str);
end
[rows cols] = size(rk.A);
x = repmat('%5.16E\t',1,(cols-1));
fprintf(writeFid,[x,'%5.16E\n'],rk.A');

str = 'b';
fprintf(writeFid,'\n%s\r\n',str);
[rows cols] = size(rk.b');
x = repmat('%5.16E\t',1,(cols-1));
fprintf(writeFid,[x,'%5.16E\n'],rk.b');

str = 'c^T';
fprintf(writeFid,'\n%s\r\n',str);
[rows cols] = size(rk.c');
x = repmat('%5.16E\t',1,(cols-1));
fprintf(writeFid,[x,'%5.16E\n\n'],rk.c');

if (class(1:2)=='2S' | class(1:2)=='3S')
    str = '==============================================================';
    fprintf(writeFid,'\n%s\r\n\n',str);
else
    % Write to file coefficients for low-storage formulation  
    str = 'alpha';
    fprintf(writeFid,'\n%s\r\n',str);
    [rows cols] = size(rk.alpha');
    x = repmat('%5.16E\t',1,(cols-1));
    fprintf(writeFid,[x,'%5.16E\n'],rk.alpha');
    
    str = 'beta';
    fprintf(writeFid,'\n%s\r\n',str);
    [rows cols] = size(rk.beta');
    x = repmat('%5.16E\t',1,(cols-1));
    fprintf(writeFid,[x,'%5.16E\n'],rk.beta');
    
    str = 'gamma1';
    fprintf(writeFid,'\n%s\r\n',str);
    [rows cols] = size(rk.gamma1');
    x = repmat('%5.16E\t',1,(cols-1));
    fprintf(writeFid,[x,'%5.16E\n'],rk.gamma1');
    
    str = 'gamma2';
    fprintf(writeFid,'\n%s\r\n',str);
    [rows cols] = size(rk.gamma2');
    x = repmat('%5.16E\t',1,(cols-1));
    fprintf(writeFid,[x,'%5.16E\n'],rk.gamma2');
    
    str = 'gamma3';
    fprintf(writeFid,'\n%s\r\n',str);
    [rows cols] = size(rk.gamma3');
    x = repmat('%5.16E\t',1,(cols-1));
    fprintf(writeFid,[x,'%5.16E\n'],rk.gamma3');
    
    str = 'delta';
    fprintf(writeFid,'\n%s\r\n',str);
    [rows cols] = size(rk.delta');
    x = repmat('%5.16E\t',1,(cols-1));
    fprintf(writeFid,[x,'%5.16E\n'],rk.delta');
    
    str = '==============================================================';
    fprintf(writeFid,'\n%s\r\n\n',str);
end

wf= 1;

