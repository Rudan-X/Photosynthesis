function [c,ceq] = mycon(x,kms)
c=[];
ceq=[];

ceq(1)=x(1)*x(3)*x(5)*x(7)*x(9)/(x(2)*x(4)*x(6)*x(8)*x(10))-1;

ceq(2)=x(5)*x(7)*x(9)/x(1)/(x(6)*x(9)+x(7)*x(9)+x(5)*x(9)+x(5)*x(7))-kms(1);

if ~isnan(kms(2))
    ceq(3)=x(9)*(x(4)*x(6)+x(4)*x(7)+x(5)*x(7))/x(3)/(x(6)*x(9)+x(7)*x(9)+x(5)*x(9)+x(5)*x(7))-kms(2);
else 
    ceq(3)=0;
end

if ~isnan(kms(3))
    ceq(4)=x(2)*(x(4)*x(6)+x(4)*x(7)+x(5)*x(7))/x(8)/(x(4)*x(6)+x(2)*x(6)+x(2)*x(4)+x(2)*x(5))-kms(3);
end

if ~isnan(kms(4))
    ceq(5)=x(2)*x(4)*x(6)/x(10)/(x(4)*x(6)++x(2)*x(6)+x(2)*x(4)+x(2)*x(5))-kms(4);
end