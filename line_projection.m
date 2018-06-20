%%% pneumococci diplo analysis

clear all
close all
%%% select image 

pix = 10; % pixel size in nm
N_scale = 1; % rescaling of image
[filename, pathname] = uigetfile({'*.tif';'*.png';'*.asc';'*.dat',},'Pick a MATLAB code file');
if isequal(filename,0) || isequal(pathname,0)
   disp('User pressed cancel')
   break
else
   disp(['User selected ', fullfile(pathname, filename)])
end
file = fullfile(pathname, filename);
disp(filename)

image = importdata(file);
image = imresize(image,N_scale);

% show the image
figure, imshow(image)
pointX = []; pointY = [];


%%% select bacterias to analyze manually 
i=1; dr=5; OK = true;
while true
    try
        [x,y] = ginput(1); % here get the coordinates for the click
        j=1;
        while j<i % checking if the click is close to some previous click
            if (x-pointX(j))^2+(y-pointY(j))^2 < dr^2
                pointX(j) = []; pointY(j) = []; % if so the point is removed
                i=i-1; % the list is reduced by one point
                OK = false; % the present click is not OK
            end
            j=j+1;
        end
        if OK % if the click is not to close any previous click it is OK
            pointX(i)=x; % point is added to the list of points
            pointY(i)=y;
            i=i+1;
        end
        OK = true; % standard assumptions is that a click is OK 
        imshow(image)
        hold on, plot(pointX,pointY,'.g','MarkerSize',18), hold off
        % Draw lines between points
        k=1; xp={}; yp={};
        while 2*k<=i-1
            xp{k} = [pointX(2*k), pointX(2*k-1)];
            yp{k} = [pointY(2*k),pointY(2*k-1)];
            hold on, plot(xp{k},yp{k},'-g')
            hold off
            k=k+1;
        end
            
    catch
        break
    end
end

figure, imshow(image)
hold on, plot(pointX,pointY,'.g','MarkerSize',18), hold off

for l=1:size(xp,2)
    hold on, plot(xp{l},yp{l},'-g'), hold off
end


%%% Separate image into two different channels (green and red), original
%%% image must have one channel as pure green, i.e 'colormap = [0 0:255 0]'
%{
[X, map] = rgb2ind(image,256);
colormap_purple = sum(map(:,3));
if colormap_purple ~= 0
    GreenIm = image(:,:,2);
    RedIm = image(:,:,1)+image(:,:,3);
else
    GreenIm = image(:,:,2);
    RedIm = image(:,:,1);
end
%}
tic
%%% convert images to double precision
GreenIm = image;
RedIm = image;
GreenIm = im2double(GreenIm);
GreenIm = GreenIm./max(max(GreenIm(:)));
RedIm = im2double(RedIm);

image = rgb2gray(image); % convert to grayscale image
image = im2double(image); % convert image to double precision
lineI={}; accI = {};
for n=1:size(xp,2)
    % start & end-points of lines
    x1 = xp{n}(1); x2 = xp{n}(2);
    y1 = yp{n}(1); y2 = yp{n}(2);
    dx = x1-x2; dy = y1-y2;
    ortI=[];
    if abs(dx)>abs(dy) % check if lines is more parallell with x-axis
        disp(['x-parallel nr. :',num2str(n)])
        if x1>x2
            xx = x2:x1;
        else
            xx = x1:x2;
        end
        % calculate the equation for line y=kx+m
        k = (y1-y2)/(x1-x2);
        m = y1-k*x1;
        % calculate the line itself
        yy = k*xx+m;
        % make sure the points on the line are integers
        %xx = round(xx); 
        %yy = round(yy);
        % find orthogonal line ort = qx+r
        q = -1/k;
        r = yy-q*xx;
        xxO=[]; yyO=[];
        green_UP_projection=zeros(1,round(0.5*length(xx))); red_UP_projection=zeros(1,round(0.5*length(xx)));
        green_DOWN_projection=zeros(1,round(0.5*length(xx))); red_DOWN_projection=zeros(1,round(0.5*length(xx)));
        for i=1:length(xx)
            xxO(i,:) = (xx(i)-8*pix*N_scale*cos(atan(q)):1:xx(i)+8*pix*N_scale*cos(atan(q)));
            yyO(i,:) = q*xxO(i,:)+r(i);
            xxO = round(xxO);
            yyO = round(yyO);
            hold on, plot(xxO(i,:),yyO(i,:),'.y')
        end
        %%% draw line in nm
        line_length = 10*sqrt((xx(end)-xx(1)).^2+(yy(end)-yy(1)).^2);
        x_line = linspace(0,line_length,size(xx,2))';
            %%% Calculate cross correlation along orthogonal lines
            for jj=1:size(xxO,1)
                green_UP_not_counted = true;
                red_UP_not_counted = true;
                green_DOWN_not_counted = true;
                red_DOWN_not_counted = true;
                for kk=1:round(0.5*size(xxO,2))
                    % if-part to avoid end up outside image boundarys
                    if yyO(jj,kk)>0 && yyO(jj,kk)<size(GreenIm,1)...
                            && xxO(jj,kk)>0 && xxO(jj,kk)<size(GreenIm,2)
                        green_patch_UP(jj,kk) = GreenIm(yyO(jj,kk),xxO(jj,kk));
                        if GreenIm(yyO(jj,kk),xxO(jj,kk))>0.2 && green_UP_not_counted
                            green_UP_projection(jj) = 1; 
                            green_UP_not_counted = false;
                        end
                        if RedIm(yyO(jj,kk),xxO(jj,kk))~=0 && red_UP_not_counted
                            red_DOWN_projection(jj) = 1;
                            red_UP_not_counted = false; 
                        end
                        if xxO(jj,end-kk)<size(GreenIm,2) && yyO(jj,end-kk)<size(GreenIm,1)...
                                && xxO(jj,end-kk)>0 && yyO(jj,end-kk)>0
                            green_patch_DOWN(jj,kk) = GreenIm(yyO(jj,end-kk),xxO(jj,end-kk)); 
                            if GreenIm(yyO(jj,end-kk),xxO(jj,end-kk))>0.2 && green_DOWN_not_counted 
                                green_DOWN_projection(jj) = 1; 
                                green_DOWN_not_counted = false;
                            end
                        end
                        if xxO(jj,end-kk)<size(RedIm,2) && yyO(jj,end-kk)<size(RedIm,1)...
                                && xxO(jj,end-kk)>0 && yyO(jj,end-kk)>0
                            if RedIm(yyO(jj,end-kk),xxO(jj,end-kk))~=0 && red_DOWN_not_counted 
                                red_DOWN_projection(jj) = 1;
                                red_DOWN_not_counted = false; 
                            end
                        end
                    end
                end
            end


   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else % else line is more parallell with y-axis
        disp(['y-parallel nr. :',num2str(n)])
        if y1>y2
            yy = y2:y1;
        else
            yy = y1:y2;
        end
        % calculate the equation for line
        k = (x1-x2)/(y1-y2);
        m = x1-k*y1;
        % calculate the line itself
        xx = k*yy+m;
        % make sure the points on the line are integers
        xx = round(xx); 
        yy = round(yy);
        % find orthogonal lines
        q = -1/k;
        r = xx-q*yy;
        xxO=[]; yyO=[];
        green_UP_projection=zeros(1,round(0.5*length(xx))); red_UP_projection=zeros(1,round(0.5*length(xx)));
        green_DOWN_projection=zeros(1,round(0.5*length(xx))); red_DOWN_projection=zeros(1,round(0.5*length(xx)));
        for i=1:length(xx)
            yyO(i,:) = yy(i)-8*pix*N_scale*cos(atan(q)):1:yy(i)+8*pix*N_scale*cos(atan(q));
            xxO(i,:) = q*yyO(i,:)+r(i);
            yyO = round(yyO);
            xxO = round(xxO);
            hold on, plot(xxO(i,:),yyO(i,:),'.y')
        end
        %%% draw line in nm
        line_length = 10*sqrt((xx(end)-xx(1)).^2+(yy(end)-yy(1)).^2);
        x_line = linspace(0,line_length,size(xx,2))';
           %%% Calculate cross correlation along orthogonal lines
           for jj=1:size(yyO,1)
                green_UP_not_counted = true;
                red_UP_not_counted = true;
                green_DOWN_not_counted = true;
                red_DOWN_not_counted = true;
                for kk=1:round(0.5*size(yyO,2))
                    % if-part to avoid end up outside image boundarys
                    if yyO(jj,kk)>0 && yyO(jj,kk)<size(GreenIm,1)...
                            && xxO(jj,kk)>0 && xxO(jj,kk)<size(GreenIm,2)
                        green_patch_UP(jj,kk) = GreenIm(yyO(jj,kk),xxO(jj,kk)); 
                        if GreenIm(yyO(jj,kk),xxO(jj,kk))>0.2 && green_UP_not_counted
                            green_UP_projection(jj) = 1; 
                            green_UP_not_counted = false;
                        end
                        if RedIm(yyO(jj,kk),xxO(jj,kk))~=0 && red_UP_not_counted
                            red_DOWN_projection(jj) = 1;
                            red_UP_not_counted = false; 
                        end
                        if xxO(jj,end-kk)<size(GreenIm,2) && yyO(jj,end-kk)<size(GreenIm,1)...
                                && xxO(jj,end-kk)>0&& yyO(jj,end-kk)>0
                            green_patch_DOWN(jj,kk) = GreenIm(yyO(jj,end-kk),xxO(jj,end-kk)); 
                            if GreenIm(yyO(jj,kk),xxO(jj,end-kk))>0.2 && green_DOWN_not_counted
                                green_DOWN_projection(jj) = 1; 
                                green_DOWN_not_counted = false;
                            end
                        end
                        if xxO(jj,end-kk)<size(RedIm,2) && yyO(jj,end-kk)<size(RedIm,1)...
                                && xxO(jj,end-kk)>0 && yyO(jj,end-kk)>0
                            if RedIm(yyO(jj,kk),xxO(jj,end-kk))~=0 && red_DOWN_not_counted 
                                red_DOWN_projection(jj) = 1;
                                red_DOWN_not_counted = false; 
                            end
                        end
                    end
                end
           end

        %end
        

    end
    %accI{n} = ortI;
    hold on, plot(xx,yy,'.b'),hold off
    projection = sum(green_patch_UP')+sum(green_patch_DOWN');
    projection = projection./sum(projection(:));
    HoM2(n) = sum(projection.^2);
    HoM4 = sum(projection.^4);
    HoM6 = sum(projection.^6);
    HoM8 = sum(projection.^8);
    %figure(2)
    %plot(x_line,projection)
    green_projection_total(n) = sum(green_UP_projection)+sum(green_DOWN_projection);
    green_area_fraction(n) = (sum(green_UP_projection)+sum(green_DOWN_projection))./2/length(green_UP_projection);
    
    bacteria_length(n) = line_length;
    red_projection_total(n) = sum(red_UP_projection)+sum(red_DOWN_projection);
    red_area_fraction(n) = (sum(red_UP_projection)+sum(red_DOWN_projection))./2/length(red_UP_projection);
end


if size(green_area_fraction,1)==1
    green_area_fraction = green_area_fraction';
end

if size(bacteria_length,1)==1
    bacteria_length = bacteria_length';
end

if size(HoM2,1)==1
    HoM2 = HoM2';
end

