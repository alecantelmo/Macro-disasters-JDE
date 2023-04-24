function result=bigprod(bigmat,smallmat)
%
% © Copyright, Oren Levintal, June 13, 2016.

try
    result=bigmat*smallmat;
catch
    m=size(bigmat,1);
    n=size(smallmat,2);
    if nnz(bigmat)==0 || nnz(smallmat)==0
        result=sparse(m,n);
    else
        [i,j,s]=find(bigmat);
         clear bigmat
        
        temp=repmat(s,1,n).*smallmat(j,:); clear j s smallmat
        [sorti,indI] = sort(i,1); clear i
        Di=[(abs(diff(sorti))~=0);1];
        cumtemp=cumsum(temp(indI,:)); clear temp indI
        sumtemp=cumtemp(Di==1,:); clear cumtemp
        sumtemp=[sumtemp(1,:);diff(sumtemp)];
        Di=[1;Di(1:end-1)];
        newi=repmat(sorti(Di==1),n,1);
        newj=repmat(1:n,length(sorti(Di==1)),1); clear sorti Di
        newj=newj(:);
        result=sparse(newi,newj,sumtemp(:),m,n);
    end
end

end