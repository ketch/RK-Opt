function wf=writeFile(A,b,c,alpha,beta,gamma1,gamma2,gamma3,delta,s,p)
%function wf=writeFile(A,b,c,alpha,beta,gamma1,gamma2,gamma3,delta,s,p)
%
% 
% Write to file Butcher's coefficients and low-storage coefficients if 
% required.

outputFileName = strcat('ERK-',num2str(p),'-',num2str(s),'.txt');
writeFid = fopen(outputFileName,'w');

fprintf(writeFid, '%s\t\t %s\n', '#stage','order');
output = [s;p];
fprintf(writeFid, '%u\t \t\t%u\n\n',output);

% Butcher coefficients
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

if isempty(alpha)
    str = '==============================================================';
    fprintf(writeFid,'\n%s\r\n\n',str);
else
    % Write to file coefficients for low-storage formulation  
    str = 'alpha';
    fprintf(writeFid,'\n%s\r\n',str);
    [rows cols] = size(alpha');
    x = repmat('%5.16E\t',1,(cols-1));
    fprintf(writeFid,[x,'%5.16E\n'],alpha');
    
    str = 'beta';
    fprintf(writeFid,'\n%s\r\n',str);
    [rows cols] = size(beta');
    x = repmat('%5.16E\t',1,(cols-1));
    fprintf(writeFid,[x,'%5.16E\n'],beta');
    
    str = 'gamma1';
    fprintf(writeFid,'\n%s\r\n',str);
    [rows cols] = size(gamma1');
    x = repmat('%5.16E\t',1,(cols-1));
    fprintf(writeFid,[x,'%5.16E\n'],gamma1');
    
    str = 'gamma2';
    fprintf(writeFid,'\n%s\r\n',str);
    [rows cols] = size(gamma2');
    x = repmat('%5.16E\t',1,(cols-1));
    fprintf(writeFid,[x,'%5.16E\n'],gamma2');
    
    str = 'gamma3';
    fprintf(writeFid,'\n%s\r\n',str);
    [rows cols] = size(gamma3');
    x = repmat('%5.16E\t',1,(cols-1));
    fprintf(writeFid,[x,'%5.16E\n'],gamma3');
    
    str = 'delta';
    fprintf(writeFid,'\n%s\r\n',str);
    [rows cols] = size(delta');
    x = repmat('%5.16E\t',1,(cols-1));
    fprintf(writeFid,[x,'%5.16E\n'],delta');
    
    str = '==============================================================';
    fprintf(writeFid,'\n%s\r\n\n',str);
end
wf= 1;

