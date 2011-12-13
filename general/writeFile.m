function wf=writeFile(rk)
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

values = struct2cell(rk);
names  = fieldnames(rk);
for i=1:length(values)
    writeField(writeFid,names{i},values{i});
end

str = '==============================================================';
fprintf(writeFid,'\n%s\r\n\n',str);

wf= 1;





