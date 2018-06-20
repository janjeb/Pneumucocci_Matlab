function [binIm] = AdaptiveThreshold(Im)
    [Lx, Ly] = size(Im);
    binIm = zeros(Lx,Ly);
    Lb = 25; % size of box in which the image is divided for thresholding
    x = fix(Lx/Lb); y = fix(Ly/Lb);
    for i=0:x-1;
        for j=0:y-1;
            sub = Im(Lb*i+1:Lb*(i+1),Lb*j+1:Lb*(j+1)); 
            level = graythresh(sub);
            subBin = im2bw(sub,level);
            binIm(Lb*i+1:Lb*(i+1),Lb*j+1:Lb*(j+1)) = subBin;
        end
    end
end