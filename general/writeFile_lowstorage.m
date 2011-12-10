function wfls=writeFile_lowstorage(rk.alpha,rk.beta,rk.gamma1,rk.gamma2,rk.gamma3,rk.delta,class,writeFid)
%function wfls=writeFile_lowstorage(rk.alpha,rk.beta,rk.gamma1,rk.gamma2,rk.gamma3,rk.delta,class,writeFid)
%
% 
% Write to file coefficients for low storage formulation

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

wfls = 1;
