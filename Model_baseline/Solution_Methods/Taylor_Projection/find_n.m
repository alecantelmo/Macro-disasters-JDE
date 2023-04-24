function n=find_n(pi_,u)
% Given u=pi_(v,u), how many compositions of type pi_(v,pi_(v,u)) are needed to get a
% function of v only.
%
% © Copyright, Oren Levintal, June 13, 2016.

if isempty(pi_)
    n=0;
else
    check=0;
    subspi_=pi_;
    n=0;
    while check==0
        n=n+1;
        newsubspi_=subs(subspi_,u,pi_);
%         check=nnz(1-logical(jacobian(subspi_,u)==0));
        check=isequal(newsubspi_,subspi_);
        subspi_=newsubspi_;
        
    end
end

end