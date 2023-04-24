function disp_deriv( fun_string,i,fun_sym , varargin)
%disp_deriv breaks derivativs into several lines if the length of derivative is
%higher than 1000 characters.
%
% © Copyright, Oren Levintal, June 13, 2016.

if isempty(varargin)
    texp=char(fun_sym(i));
else
    texp=char(fun_sym(varargin{1}));
end
if length(texp)<1000
    disp([fun_string '(' num2str(i)  ')=' texp ';']);
    texp='';
else
    itexp=1000;
    while strcmp(texp(itexp),'+')+strcmp(texp(itexp),'-')+strcmp(texp(itexp),'*')+strcmp(texp(itexp),'/')==0 & itexp<length(texp)
        itexp=itexp+1;
    end
    if itexp==length(texp)
        disp([fun_string '(' num2str(i)  ')=' texp(1:itexp) ';' ]);
        texp='';
    else
        disp([fun_string '(' num2str(i)  ')=' texp(1:itexp) '...' ]);
        texp=texp(itexp+1:end);
    end
end

while length(texp)>0
    itexp=min(1000,length(texp));
    while strcmp(texp(itexp),'+')+strcmp(texp(itexp),'-')+strcmp(texp(itexp),'*')+strcmp(texp(itexp),'/')==0 & itexp<length(texp)
        itexp=itexp+1;
    end
    if itexp==length(texp)
        disp([texp ';']);
        texp='';
    else
        disp([texp(1:itexp) '...' ]);
        texp=texp(itexp+1:end);
    end
end
R=0;

end

