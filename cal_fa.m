% l1 = 1.2;
% l2 = 0.4;
% l3 = 0.4;

function fa = cal_fa(l1,l2,l3)
%[0.00084860441 , 0.0011277209 ,0.00063587740]
% l1 =     1.6136  
%     l2 =      0.9078
%     l3 = 0.6224
md = (l1 + l2 + l3) / 3;

fa = sqrt(3/2) * sqrt((l1 - md).^2 + (l2 - md).^2 + (l3 - md).^2) ./ (sqrt(l1.^2 + l2.^2 + l3.^2));