function [ y ] = intarray( x )
%Convert to int32 or int64, depending on the operating system.
%
% © Copyright, Oren Levintal, June 13, 2016.

archstr = computer('arch');
if strcmp(archstr(end-1:end),'32')
    y=int32(x);
elseif strcmp(archstr(end-1:end),'64')
    y=int64(x);
else
    error('cannot find if MATLAB runs on 32-bit or 64-bit operating system.')
end
end

