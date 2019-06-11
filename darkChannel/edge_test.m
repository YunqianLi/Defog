function mark = edge_test(I)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[dimr dimc c]=size(I); 
if(c==3)
    Im=rgb2gray(I);
Edge_Dg=edge(Im,'prewitt',0.05);
imshow(Edge_Dg);
mark=sum(sum(Edge_Dg,1),2)/(dimr*dimc);


end

