function gen_fun(f,symparams,v,fname,varargin)
%
% © Copyright, Oren Levintal, June 13, 2016.

rowformat=0; % function returns a column vector (number of columns>1 in vectorized form)
if ~isempty(varargin)
    if strcmp(varargin{1},'row')
        rowformat=1; % function returns a row vector (number of rows>1 in vectorized form)
    end
end

f=f(:);
if numel(f)>1
    nnzf=1-logical(f==0);
    index=find(nnzf);
else
    index=1;
end
fid = fopen([fname '_fun.m'], 'w');
fprintf(fid,'%s\n', ['function function_f=' fname '_fun(variables_v,parameters)']);
for i=1:length(symparams)
    if strcmp(char(symparams(i)),'function_f') || strcmp(char(symparams(i)),'variables_v') || strcmp(char(symparams(i)),'parameters')
        error([char(symparams(i)) 'is reserved. Change variable/parameter name.']);
    end
    fprintf(fid,'%s\n', [char(symparams(i)) '=parameters(' num2str(i) ');']);
end
for i=1:length(v)
    if strcmp(char(v(i)),'function_f') || strcmp(char(v(i)),'variables_v') || strcmp(char(v(i)),'parameters')
        error([char(v(i)) 'is reserved. Change variable/parameter name.']);
    end
    fprintf(fid,'%s\n', [char(v(i)) '=variables_v(' num2str(i) ',:);']);
end
if isempty(f)
    fprintf(fid,'%s\n', ['function_f=[];']);
else
    if rowformat==0
        fprintf(fid,'%s\n', ['function_f=zeros(' num2str(numel(f)) ',size(variables_v,2));']);
    elseif rowformat==1
        fprintf(fid,'%s\n', ['function_f=zeros(size(variables_v,2),' num2str(numel(f)) ');']);
    end
    for i=index(:)'
        if rowformat==0
            disp_fun('function_f',f(i),i,fid);
        elseif rowformat==1
            disp_fun_row('function_f',f(i),i,fid);
        end
    end
end
fclose(fid);
