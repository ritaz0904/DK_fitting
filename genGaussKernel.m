function f = genGaussKernel(vox,FWHM)

pixFWHM = FWHM./vox;
s = pixFWHM./(2*sqrt(2*log(2)));
mxSize = ceil(s);

scale = round(3);
f = zeros(2*scale*mxSize(1)+1,2*scale*mxSize(2)+1,2*scale*mxSize(3)+1);
for x = -scale*mxSize(1):scale*mxSize(1)
    for y = -scale*mxSize(2):scale*mxSize(2)
        for z = -scale*mxSize(3):scale*mxSize(3)
            A = 1/sqrt(2*pi*s(1)^2)*1/sqrt(2*pi*s(2)^2)*1/sqrt(2*pi*s(3)^2);
            temp = A*exp(-(x^2/(2*s(1)^2)+y^2/(2*s(2)^2)+z^2/(2*s(3)^2)));
            f(x+scale*mxSize(1)+1,y+scale*mxSize(2)+1,z+scale*mxSize(3)+1) = temp;
        end
    end
end