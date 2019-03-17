
function ROI=findROI(Image,Radius,Center)

[ImgW,ImgH] = size(Image);

x = linspace(0, 2*pi, 50);

mask = poly2mask(Radius*cos(x)+Center(1), Radius*sin(x)+Center(2), ImgW, ImgH);

ROI = immultiply(Image,mask);

