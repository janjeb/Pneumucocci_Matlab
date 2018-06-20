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
file = fullfile(pathname, filename)

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
[X, map] = rgb2ind(image,256);
colormap_purple = sum(map(:,3));
if colormap_purple ~= 0
    GreenIm = image(:,:,2);
    RedIm = image(:,:,1)+image(:,:,3);
else
    GreenIm = image(:,:,2);
    RedIm = image(:,:,1);
end

tic
%%% convert images to double precision
GreenIm = im2double(GreenIm); 
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
        xx = round(xx); 
        yy = round(yy);
        
        green_profile = zeros(size(xx,2),1);
        red_profile = zeros(size(xx,2),1);
        for i=1:size(xx,2)
            green_profile(i) = GreenIm(yy(i),xx(i));
            red_profile(i) = RedIm(yy(i),xx(i));
        end
        line_length = 10*sqrt((xx(end)-xx(1)).^2+(yy(end)-yy(1)).^2);
        RG_dist(n) = line_length;
        x_line = linspace(0,line_length,size(xx,2))';

        fittype = 'a*exp(-((x-b)/c)^2)';
        options = fitoptions(fittype);
        options.StartPoint = [1,0.5*line_length,100];
        options.Lower = [0,0,10];
        options.Upper = [10,line_length,500];
        green_fit = fit(x_line,green_profile,fittype,options);
        red_fit = fit(x_line,red_profile,fittype,options);
        
        RG_gap(n) = abs(green_fit.b-red_fit.b);
        figure
        plot(green_fit,x_line,green_profile,'g')
        hold on
        plot(red_fit, x_line, red_profile,'r')
        hold off
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
        
        green_profile = zeros(size(xx,2),1);
        red_profile = zeros(size(xx,2),1);
        for i=1:size(xx,2)
            green_profile(i) = GreenIm(yy(i),xx(i));
            red_profile(i) = RedIm(yy(i),xx(i));
        end
        line_length = 10*sqrt((xx(end)-xx(1)).^2+(yy(end)-yy(1)).^2);
        RG_dist(n) = line_length;
        x_line = linspace(0,line_length,size(xx,2))';
        
        fittype = 'a*exp(-((x-b)/c)^2)';
        options = fitoptions(fittype);
        options.StartPoint = [1,0.5*line_length,40];
        options.Lower = [0,0,10];
        options.Upper = [10,line_length,500];
        green_fit = fit(x_line,green_profile,fittype,options);
        red_fit = fit(x_line,red_profile,fittype,options);
        
        RG_gap(n) = abs(green_fit.b-red_fit.b);
        figure
        plot(green_fit, x_line, green_profile,'g')
        hold on
        plot(red_fit, x_line, red_profile,'r')
        hold off
    end
    %accI{n} = ortI;
    %hold on, plot(xx,yy,'.b'),hold off

end
RG_dist = RG_dist';
RG_gap = RG_gap';








