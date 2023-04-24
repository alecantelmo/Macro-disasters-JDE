function [ newf ] = subsf( f,auxvars,auxfuns )
%Substitute out auxiliary variables.
%
% © Copyright, Oren Levintal, June 13, 2016.

f_temp=[]; newf=f;
while ~isequal(f_temp,newf)
    f_temp=newf;
    newf=subs(newf,auxvars,auxfuns);
end

end

