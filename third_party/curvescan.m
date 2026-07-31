function curvescan(datadir,graphimagefile,graphdatafile)
% curvescan(graphimagefile,graphdatafile)
% This function reads in the image file of a graph and allows the user to
% extract data from a curve on the graph.  The user sets the coordinate
% transformation from image to data coordinates with a series of mouse
% clicks on the graph image.  The curve data are calculated using cubic
% spline interpolation of user mouse clicks on the curve in the graph
% image.
%
% Jeremy Teichman, Institute for Defense Analyses, August 26, 2009

% Read in the image
if nargin == 0 | nargin == 1
    if nargin == 1  cd(datadir); end
	[FileName,PathName] = uigetfile('*.*','Load curve image input file');
	if FileName ~= 0
		graphimage = imread(strcat(PathName,FileName));
	else
		display('No curve image');
		return;
    end
else
	graphimage = imread(graphimagefile);
end

interp_points = 100; % Number of points to use in the interpolation

program_name = 'curvescan'; % For title bar of dialog boxes

figure(1);
clf;
image(graphimage); % Plot the graph image

% Collect information for coordinate transform from image to data
h = msgbox('Click origin and head of x-axis of plot',program_name);
uiwait(h);
[xbound,ybound] = ginput(2);

% In an image 0,0 is the upper right
[a,LL] = min(xbound-ybound); % Lower left
[a,LR] = max(xbound+ybound); % Lower right

Rot = [[(xbound(LR)-xbound(LL)),-(ybound(LR)-ybound(LL))]/...
	sqrt((xbound(LR)-xbound(LL))^2+(ybound(LR)-ybound(LL))^2)];

Rot = [Rot ; [-Rot(2) Rot(1)]];

h = msgbox('Click anchor point or press Enter to use origin',program_name);
uiwait(h);
[x0prime,y0prime] = ginput(1);
if isempty(x0prime)
	x0prime = xbound(LL);
	y0prime = ybound(LL);
end

origin = inputdlg({'X-coordinate of anchor point','Y-coordinate of anchor point'},program_name);
x0 = str2double(origin(1));
y0 = str2double(origin(2));

h = msgbox('Click x-axis interval',program_name);
uiwait(h);
[x_xscale y_xscale] = ginput(2);

xinterval = str2double(inputdlg('X-axis interval',program_name));
xscale = xinterval/sqrt((x_xscale(2)-x_xscale(1))^2+(y_xscale(2)-y_xscale(1))^2);

h = msgbox('Click y-axis interval',program_name);
uiwait(h);
[x_yscale y_yscale] = ginput(2);

yinterval = str2double(inputdlg('Y-axis interval',program_name));
yscale = yinterval/sqrt((x_yscale(2)-x_yscale(1))^2+(y_yscale(2)-y_yscale(1))^2);
ScaleRot = diag([xscale,yscale])*Rot;


xmin = ScaleRot*[xbound(LL)-x0prime , -(ybound(LL)-y0prime)]';
xmin = xmin(1)+x0;
xmax = ScaleRot*[xbound(LR)-x0prime , -(ybound(LR)-y0prime)]';
xmax = xmax(1)+x0;
xx = [xmin:(xmax-xmin)/interp_points:xmax];

% Transform image into data coordinates for comparison to data curve
intensity = sum(graphimage,3);
imdim = size(intensity);
[xdim,ydim] = meshgrid([1:imdim(2)],[1:imdim(1)]);
xdim2prime = xdim(intensity<50);
ydim2prime = ydim(intensity<50);
xdim2 = ScaleRot*[xdim2prime-x0prime , -(ydim2prime-y0prime)]';
xdim2(1,:) = xdim2(1,:) + x0;
xdim2(2,:) = xdim2(2,:) + y0;

% Loop until user is satisfied with interpolated data curve
okflag = 0;
while ~okflag
	figure(1);
	h = msgbox('Click points on curve, press Enter when done',program_name);
	uiwait(h);
	[xprime,yprime] = ginput();

	x = ScaleRot*[xprime-x0prime , -(yprime-y0prime)]';
	x(1,:) = x(1,:) + x0;
	x(2,:) = x(2,:) + y0;

	yy = interp1(x(1,:),x(2,:),xx,'spline');

	% Plot graph image (blue) and data curve (red)
	figure(2);
	clf;
	plot(xdim2(1,:),xdim2(2,:),'b.','MarkerSize',1);
	hold on;
	plot(xx,yy,'r');

	okstr = questdlg('Keep this curve?',program_name,'Yes','No','Yes');
	okflag = strcmp(okstr,'Yes') || strcmp(okstr,'');
end
	
	graphsurf = yy;
	figure(1);

% Send output to file (or to command window if none is chosen
if nargin > 1
	outfile = fopen(graphdatafile,'wt');
else
	[FileName,PathName] = uiputfile('*.*','Curve data output file');
	if FileName ~= 0
		outfile = fopen(strcat(PathName,FileName),'wt');
	else
		outfile = 1; % Standard Output
	end
end

for j = 1:100
	fprintf(outfile,'%f %f\n',xx(j),graphsurf(j));
end
	
if outfile ~= 1
	fclose(outfile);
end