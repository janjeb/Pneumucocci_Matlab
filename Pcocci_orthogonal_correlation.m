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
        % find orthogonal line ort = qx+r
        q = -1/k;
        r = yy-q*xx;
        xxO=[]; yyO=[];
        green_patch=[]; red_patch=[];
        for i=1:length(xx)
            xxO(i,:) = (xx(i)-8*pix*N_scale*cos(atan(q)):0.1:xx(i)+8*pix*N_scale*cos(atan(q)));
            yyO(i,:) = q*xxO(i,:)+r(i);
            xxO = round(xxO);
            yyO = round(yyO);
            hold on, plot(xxO(i,:),yyO(i,:),'.y')
            
            %%% Calculate cross correlation along orthogonal lines
            for jj=1:size(xxO,1)
                for kk=1:size(xxO,2)
                    % values along orthogonal line j
                    %Ggreen(k) = GreenIm(xxO(j,k),yyO(j,k)); 
                    %Gred(k) = RedIm(xxO(j,k),yyO(j,k));
                    green_patch(jj,kk) = GreenIm(yyO(jj,kk),xxO(jj,kk)); 
                    red_patch(jj,kk)   = RedIm(yyO(jj,kk),xxO(jj,kk));
                end
            end
            %%% Calculate cross correlation
            %{
            % cross corelation as average of red-green and green-red
            for row=1:size(red_patch,1)
                Ggreen = green_patch(row,:);
                Gred   = red_patch(row,:);
                %%% Calculate cross correlation
                end_pt = round(0.5*size(xxO,2)); % end point of ortogonal line
                % start point for correlation
                start_pt = round(0.5*size(Gred,2));
                if mod(size(Gred,2),2)==0
                    start_pt1 = start_pt+1;
                    start_pt2 = start_pt;
                else
                    start_pt1 = start_pt;
                    start_pt2 = start_pt;
                end
                Ggreen_mean = mean(Ggreen); Gred_mean = mean(Gred);
                for ii=0:end_pt-1
                    Crg_up = 0; Cgr_up = 0; Cgr_down = 0; Crg_down = 0;
                    for jj=0:end_pt-ii-1
                        % cross-correlation along symmetry line
                        Crg_up   = (Gred(start_pt1+jj+ii)-Gred_mean).*(Ggreen(jj+start_pt1)-Ggreen_mean) + Crg_up;
                        Cgr_up   = (Gred(jj+start_pt1)-Gred_mean).*(Ggreen(start_pt1+jj+ii)-Ggreen_mean) + Cgr_up;
                        Crg_down = (Gred(start_pt2-jj-ii)-Gred_mean).*(Ggreen(start_pt2-jj)-Ggreen_mean) + Crg_down;
                        Cgr_down = (Gred(start_pt2-jj)-Gred_mean).*(Ggreen(start_pt2-jj-ii)-Ggreen_mean) + Cgr_down;
                    end
                    % normalize cross correlation
                    CCrg_up(ii+1) = Crg_up/(end_pt-ii);
                    CCgr_up(ii+1) = Cgr_up/(end_pt-ii);
                    CCrg_down(ii+1) = Crg_down/(end_pt-ii);
                    CCgr_down(ii+1) = Cgr_down/(end_pt-ii);
                end
                % final cross correlation
                CC_up(row,:) = (CCrg_up + CCgr_up)./2;
                CC_down(row,:) = (CCrg_down + CCgr_down)./2;
            end
            %}
            %%% sum along orthogonal lines
            %acci=0;
            %for j=1:length(xxO)
            %    acci = image(yyO(j),xxO(j)) + acci;
               %fprintf(fid,'%4f%\r\n',image(yyO(j),xxO(j)));
        end
            %fprintf(fid,'bacteria \r\n');
            %ortI(i) = acci;
        %end
        for pp = 1:size(green_patch,1)
            for qq = 1:size(green_patch,2)
                green_pix_intensity(pp,qq) = green_patch(pp,qq);
                red_pix_intensity(pp,qq) = red_patch(pp,qq);
            end
        end
        green_intensity{n} = green_pix_intensity(:); 
        red_intensity{n}   = red_pix_intensity(:);
        G_mean_intensity(n) = mean(green_pix_intensity(:));
        G_std_intensity(n) = std(green_pix_intensity(:));
        R_mean_intensity(n) = mean(red_pix_intensity(:));
        R_std_intensity(n) = std(red_pix_intensity(:));
        %green_m = sum(green_patch)./size(green_patch,1);
        %red_m = sum(red_patch)./size(red_patch,1);
        
        %L_orto = pix/N_scale*sqrt((xxO(1,1)-xxO(1,end)).^2+(yyO(1,1)-yyO(1,end)).^2); % length of orthogonal line in pixels
        %base_line = 0:L_orto/(size(xxO,2)):L_orto;
        %base_line = base_line(1:end-1);
        %base_line = linspace(0,L_orto,size(xxO,2));
        %MPD = round(size(base_line,2)/2);
        %[pks_green, locs_green] = findpeaks(green_m,'MINPEAKDISTANCE',MPD);
        %[pks_red, locs_red] = findpeaks(red_m,'MINPEAKDISTANCE',MPD);
        % width of red (units of 10 nm)
       %red_width = base_line(locs_red(2))-base_line(locs_red(1)); 
        % width of green (units of 10 nm)
       %green_width = base_line(locs_green(2))-base_line(locs_green(1));
        % difference i width divided by two. (units of nm)
       %overlap(n) = (red_width-green_width)/2; %'overlap' in nm for bacteria n     
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
        green_patch=[]; red_patch=[];
        for i=1:length(xx)
            yyO(i,:) = yy(i)-8*pix*N_scale*cos(atan(q)):0.1:yy(i)+8*pix*N_scale*cos(atan(q));
            xxO(i,:) = q*yyO(i,:)+r(i);
            yyO = round(yyO);
            xxO = round(xxO);
            hold on, plot(xxO(i,:),yyO(i,:),'.y')
           %%% Calculate cross correlation along orthogonal lines
           for jj=1:size(xxO,1)
                for kk=1:size(xxO,2)
                    % values along orthogonal line j
                    %Ggreen(k) = GreenIm(xxO(j,k),yyO(j,k)); 
                    %Gred(k) = RedIm(xxO(j,k),yyO(j,k));
                    green_patch(jj,kk) = GreenIm(yyO(jj,kk),xxO(jj,kk)); 
                    red_patch(jj,kk)   = RedIm(yyO(jj,kk),xxO(jj,kk));
                end
           end
            %{
            for j=1:size(xxO,1)
                for k=1:size(xxO,2)
                    % values along orthogonal line j
                    Ggreen(k) = GreenIm(yyO(j,k),xxO(j,k)); 
                    Gred(k) = RedIm(yyO(j,k),xxO(j,k));
                end
            end
            %%% Calculate cross correlation
            end_pt = size(xxO,2); % end point of ortogonal line
            % cross corelation as average of red-green and green-red
            for ii=0:end_pt-1
                Crg = 0; Cgr = 0;
                for jj=1:end_pt-ii
                    % cross-correlation along symmetry line
                    Crg = (Gred(jj+ii)-mean(Gred)).*(Ggreen(jj)-mean(Ggreen)) + Crg;
                    Cgr = (Gred(jj)-mean(Gred)).*(Ggreen(jj+ii)-mean(Ggreen)) + Cgr;
                end
                % normalize cross correlation
                CCrg(ii+1) = Crg/(end_pt-ii);
                CCgr(ii+1) = Cgr/(end_pt-ii);
            end
            % final cross correlation
            CC(i,:) = (CCrg + CCgr)./2;
            %}
            %%% sum along orthogonal lines
            %acci=0;
            %for j=1:length(yyO)
            %    acci = image(yyO(j),xxO(j)) + acci;
            %    %fprintf(fid,'%4f\r\n',image(yyO(j),xxO(j)));
            %end
            %fprintf(fid,'bacteria \r\n');
            %ortI(i) = acci;
        end
        for pp = 1:size(green_patch,1)
            for qq = 1:size(green_patch,2)
                green_pix_intensity(pp,qq) = green_patch(pp,qq);
                red_pix_intensity(pp,qq) = red_patch(pp,qq);
            end
        end
        green_intensity{n} = green_pix_intensity(:); 
        red_intensity{n}   = red_pix_intensity(:);
        G_mean_intensity(n) = mean(green_pix_intensity(:));
        G_std_intensity(n) = std(green_pix_intensity(:));
        R_mean_intensity(n) = mean(red_pix_intensity(:));
        R_std_intensity(n) = std(red_pix_intensity(:));
        %green_m = sum(green_patch)/size(green_patch,1);
        %red_m = sum(red_patch)/size(red_patch,1);
        
        %L_orto = pix/N_scale*sqrt((xxO(1,1)-xxO(1,end)).^2+(yyO(1,1)-yyO(1,end)).^2); % length of orthogonal line in pixels
        %base_line = 0:L_orto/(size(yyO,2)):L_orto;
        %base_line = base_line(1:end-1);
        %base_line = linspace(0,L_orto,size(yyO,2));
        %MPD = round(size(base_line,2)/2);
        %[pks_green, locs_green] = findpeaks(green_m,'MINPEAKDISTANCE',MPD);
        %[pks_red, locs_red] = findpeaks(red_m,'MINPEAKDISTANCE',MPD);
        % width of red (units of 10 nm)
       %red_width = base_line(locs_red(2))-base_line(locs_red(1)); 
        % width of green (units of 10 nm)
       %green_width = base_line(locs_green(2))-base_line(locs_green(1));
        % difference i width divided by two. (units of nm)
       %overlap(n) = (red_width-green_width)/2; %'overlap' in nm for bacteria n 
    end
    %accI{n} = ortI;
    hold on, plot(xx,yy,'.b'),hold off
    %I = [];
    %for i=1:length(xx)
    %    I(i) = image(yy(i),xx(i));% get the intensity along the line
    %    for j=length(xxO)
    %        Isum(n,i) = image(yyO(j),xxO(j));
    %    end
    %end
    %lineI{n}=I;
    %figure, plot(I)
end


G_mean_intensity = G_mean_intensity';
G_std_intensity  = G_std_intensity';
R_mean_intensity = R_mean_intensity';
R_std_intensity  = R_std_intensity';

intensity(:,1) = G_mean_intensity;
intensity(:,2) = G_std_intensity;
intensity(:,3) = R_mean_intensity;
intensity(:,4) = R_std_intensity;

toc

%overlap = overlap';
%green_m=green_m/max(green_m);
%red_m=red_m/max(red_m);
%figure, plot(base_line, green_m,'g','linewidth',2)
%hold on, plot(base_line, red_m,'r','linewidth',2)
%hold off

%Gcc = ifft(fft(green_m).*conj(fft(red_m)));
%figure, plot(Gcc,'k')
%figure, plot(CC_up')
%CC_up_mean = mean(CC_up);
%figure, plot(CC_up_mean)
%N = 2; % trace expanding exponent
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








