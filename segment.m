%parameters input
%get the radius of the circle of this slice:Radius ,and the center coordinates:Center
%get the radius of the ASC circle of previous slice:pre_ASC_R ,and the center coordinates:pre_ASC_C
%get the ASC distance of the center coordinates between the previous and the next slice:limitASC

%data output
%get the radius and the center coordinates of the ASC circle of this slice:ASC_R,ASC_C
%The ASC distance of the center is ok between the two slices:isASC=1£»if fails£¬isASC=0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ASC_R,ASC_C,isASC]=segment(Radius,Center,pre_ASC_R,pre_ASC_C,limitASC)

isASC=0;
[m,n]=size(Radius);  %get the center coordinates
minASC=9999;
ASC_R=pre_ASC_R+1.4;
ASC_C=pre_ASC_C;
if(m==0)
else
    for i=1:m
        [isOkA,distance_ASC]=distancecheck(pre_ASC_C,Center(i,:),limitASC);
        if(isOkA==0)
            
        else
            if(minASC > distance_ASC)
                ASC_C(:)=Center(i,:);
                ASC_R=Radius(i)+1.4;
                %
                isASC=1;
            end
        end
    end
end
end