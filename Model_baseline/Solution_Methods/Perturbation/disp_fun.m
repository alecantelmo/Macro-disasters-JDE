function R=disp_fun(fun_string,fun_sym,i,fid)
%disp_fun breaks a function into several lines if the length of the function is
%larger than 1000 characters.
%
% © Copyright, Oren Levintal, June 13, 2016.

texp=char(fun_sym);

% texp=strrep(texp,'*','.*');
% texp=strrep(texp,'/','./');
% texp=strrep(texp,'^','.^');
if length(texp)<1000
    R=[fun_string '(' num2str(i)  ')=' texp ';'];
    fprintf(fid,'%s\n',R);
    texp='';
else
    itexp=1000;
    while strcmp(texp(itexp),'+')+strcmp(texp(itexp),'-')+strcmp(texp(itexp),'*')+strcmp(texp(itexp),'/')==0 && itexp<length(texp)
        itexp=itexp+1;
    end
    if itexp==length(texp)
        R=[fun_string '(' num2str(i)  ')=' texp(1:itexp) ';' ];
        fprintf(fid,'%s\n',R);
        texp='';
    else
        R=[fun_string '(' num2str(i)  ')=' texp(1:itexp) '... ' ];
        fprintf(fid,'%s\n',R);
        texp=texp(itexp+1:end);
    end
end

while isempty(texp)==0
    itexp=min(1000,length(texp));
    while strcmp(texp(itexp),'+')+strcmp(texp(itexp),'-')+strcmp(texp(itexp),'*')+strcmp(texp(itexp),'/')==0 && itexp<length(texp)
        itexp=itexp+1;
    end
    if itexp==length(texp)
        R=[' ' texp ';'];
        fprintf(fid,'%s\n',R);
        texp='';
    else
        R=[' ' texp(1:itexp) '... ' ];
        fprintf(fid,'%s\n',R);
        texp=texp(itexp+1:end);
    end
end

end

