function [isOK,distance]=distancecheck(p1,p2,limit)
distance=norm(p1-p2);
if(distance > limit)
    isOK=0;
else
    isOK=1;
end
end