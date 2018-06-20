%%% pneumococci diplo analysis

clear all
%%% select image 
[filename, pathname] = uigetfile({'*.tif';'*.png';'*.asc';'*.dat',},'Pick a MATLAB code file');
if isequal(filename,0) || isequal(pathname,0)
   disp('User pressed cancel')
   break
else
   disp(['User selected ', fullfile(pathname, filename)])
end
file = fullfile(pathname, filename);

% show the image
figure, imshow(file)
pointX = []; pointY = [];


%%% select bacterias to analyze manually 
i=1; r=5; OK = true;
while true
    try
        [x,y] = ginput(1); % here get the coordinates for the click
        j=1;
        while j<i % checking if the click is close to some previous click
            if (x-pointX(j))^2+(y-pointY(j))^2 < r^2
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
        imshow(file)
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

figure, imshow(file)
hold on, plot(pointX,pointY,'.g','MarkerSize',18), hold off

for l=1:size(xp,2)
    hold on, plot(xp{l},yp{l},'-g'), hold off
end

%%% Calculate lines
image = importdata(file); 

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
        % find orthogonal line ort = qx+r
        q = -1/k;
        r = yy-q*xx;
        for i=1:length(xx)
            xxO(i,:) = xx(i)-70*cos(atan(q)):xx(i)+70*cos(atan(q));
            yyO(i,:) = q*xxO(i,:)+r(i);
            xxO = round(xxO);
            yyO = round(yyO);
            hold on, plot(xxO(i,:),yyO(i,:),'.y')
            %%% Calculate cross correlation along orthogonal lines
            for j=1:size(xxO,1)
                for k=1:size(xxO,2)
                    Ggreen(k) = GreenIm(yyO(j,k),xxO(j,k));
                    Gred(k) = RedIm(yyO(j,k),xxO(j,k));
                end
            end
            
            %%% sum along orthogonal lines
            %acci=0;
            %for j=1:length(xxO)
            %    acci = image(yyO(j),xxO(j)) + acci;
                %fprintf(fid,'%4f%\r\n',image(yyO(j),xxO(j)));
            %end
            %fprintf(fid,'bacteria \r\n');
            %ortI(i) = acci;
        end
    else % else line is more parallell with y-axis
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
        for i=1:length(xx)
            yyO(i,:) = yy(i)-70*cos(atan(q)):yy(i)+70*cos(atan(q));
            xxO(i,:) = q*yyO(i,:)+r(i);
            yyO = round(yyO);
            xxO = round(xxO);
            hold on, plot(xxO(i,:),yyO(i,:),'.y')
            %%% Calculate cross correlation along orthogonal lines
            for j=1:size(xxO,1)
                for k=1:size(xxO,2)
                    Ggreen(k) = GreenIm(yyO(j,k),xxO(j,k));
                    Gred(k) = RedIm(yyO(j,k),xxO(j,k));
                end
            end
            %%% sum along orthogonal lines
            %acci=0;
            %for j=1:length(yyO)
            %    acci = image(yyO(j),xxO(j)) + acci;
            %    %fprintf(fid,'%4f\r\n',image(yyO(j),xxO(j)));
            %end
            %fprintf(fid,'bacteria \r\n');
            %ortI(i) = acci;
        end
    end
    %accI{n} = ortI;
    hold on, plot(xx,yy,'.b'),hold off
    I = [];
    for i=1:length(xx)
        I(i) = image(yy(i),xx(i));% get the intensity along the line
        %for j=length(xxO)
        %    Isum(n,i) = image(yyO(j),xxO(j));
        %end
    end
    lineI{n}=I;
    %figure, plot(I)
end


N = 2; % trace expanding exponent
%expandedLineI = {}; expandedAccI = {}; CoM = [];
%for i=1:n
%    expandedLineI{i} = expand(lineI{i},N);
%    [expandedAccI{i},CoM(i)] = expand(accI{i},N);
%end
 

%hold on, plot(yyO,xxO,'.y')
%for n=1:size(xp,2)
%    figure, plot(expandedLineI{n})
%    figure, plot(expandedAccI{n})
%    hold on, line([CoM(n) CoM(n)], [min(expandedAccI{n}) max(expandedAccI{n})]), hold off
%end

%{
%%% save all the data to txt-file
fid = fopen([pathname,filename,'.txt'],'w');
%head = ['colour: ',color,', num granules = ',num2str(N)];

for n=1:size(xp,2)
    L =length(expandedAccI{n})*0.0025; 
    fprintf(fid,'length = %4f\r\n', L);
    fprintf(fid,'Center of mass = %4f\r\n', CoM(n));
    fprintf(fid,'%6s\r\n',['Intensity trace bacteria ',num2str(n)]);
    fprintf(fid,'%4f\r\n',expandedAccI{n});
    fprintf(fid,'\r\n');
end
%fprintf(fid,head);
%fprintf(fid,'\r\n \r\n');
%fprintf(fid,'%6s %12s\r\n','X','Y');
%fprintf(fid,'%6i %12i\r\n',pos);
fclose(fid);
%}








