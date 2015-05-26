%===================================================================
%FUNCTION: GRIDIMAGE_INTERPOLATE
%===================================================================
%   
%   [X,Y,DATA_OUT]=GRIDIMAGE_INTERPOLATE(DATA_IN,STEP_INTERP,INTERP_METHOD,OPTION)
%       Generates an interpolated image from an input image.  input
%       dimensions are m x n, output is a sucessive averaging (iterations
%       specified by STEP_INTERP) that yields an m* x n* array. Uses
%       INTERP2 function for 2D interpolation.

%------------------------------------------------------------------------------------------------------------------
%
%   INPUT
%       data = m x n matrix with (:,1) as x coodrinate, (:,2) as y and n-2 variables to image
%       step_interp = number of steps to add to x and y (must be the same)
%       interp_method = option for interpolation, where
%                       1 = nearest neighbor
%                       2 = 2D linear (default);
%                       3 = spline
%                       4 = cubic interpolation
%
%                       option = turn on/off interpolation (1=on is default)
%
%   OUTPUT
%       data_out  =  matrix of a x b for each of the variables which is the
%                           format used for imaging
%       dx,dy is step size
%altereed by DMR
%       data_mask = matrix where 0 indicated there was no data at this
%       location read from file
%
function [x,y,data_out,dx,dy,data_mask]=gridimage_interpolate(data_in,step_interp,interp_method,option,cust_bound,xbound,ybound)

if nargin<4; option=1; end
if nargin <3;interp_method=2;end;

%declarations------------------------------------------------------------------------------------------------------------------
%get new grid size axes from interpolation, use interp2 function for
%data averaging later in code
%temp=[];temp1=[];temp3=[];
%x_new=[];y_new=[]; 
data_temp=[];data_temp2=[];data_out=[];


%data reassignments--------------------------------------------------------------------------------------------------------
[npos nvar] = size(data_in);
nvar=nvar-2;    %loose distance columns
x_old = data_in(:,1);
y_old = data_in(:,2);
sw=0;

if length(x_old)<8 || length(y_old)<8
    disp('WARNING: Switching interpolation to nearest for first round because there are too few points');
    sw=1;
end


if interp_method==1;imethod='nearest';
elseif interp_method==4;imethod='cubic';
elseif interp_method==3;imethod='spline';
else imethod='linear';
end   %end checking interp options

%interpolate data matrix--------------------------------------------------------------------------------------------------------

if cust_bound %change the bounds of the data matrix, assume 10 mm spacing
    xi=[xbound(1):10:xbound(2)];
    yi=[ybound(1):10:ybound(2)];
else
    %generate new distances
    xi=unique(data_in(:,1));    %all elements in x/y axes but sorted
    yi=sortrows(unique(data_in(:,2)),-1);
end

%original matrix lengths
%nx0 = length(xi);
%ny0 = length(yi);

%altered by DMR
nx0 = (max(xi)-min(xi))/10 + 1; %assumes 10 mm spacing but will not fill in missing elements with data (will mask)
ny0 = (max(yi)-min(yi))/10 + 1;


%check spacing
if length(xi)>1
dt_x = getstep(xi);
else
    dt_x=10;
end

dt_y = getstep(yi);
if length(dt_x) > 1;disp('WARNING: X Axis not evenly spaced points');end;
if length(dt_y) > 1;disp('WARNING: Y Axis not evenly spaced points');end;

%allocate new data matricies
%data_temp = zeros(ny0,nx0);
%data_out = zeros((2^step_interp)*(ny0-1)+1,   (2^step_interp)*(nx0-1)+1,nvar);
data_temp = ones(nx0,ny0);
data_temp=data_temp.*9999;   %use nan so we know if a position has not been filled
data_out = zeros((2^step_interp)*(nx0-1)+1,   (2^step_interp)*(ny0-1)+1,nvar);
%must remake xi and yi inorder for the function ifind to work correctly
xi=[min(xi):10:max(xi)]';    %also assumes 10 mm spacing
yi=[min(yi):10:max(yi)]';


%loop to get data
for l=1:nvar
    data_temp = ones(nx0,ny0);
    data_temp=data_temp.*99999;
    
    %find the location and correctly assign in new matrix
    for i=1:npos
        x0 = data_in(i,1);  x0_idx = ifind(xi,x0);
        y0 = data_in(i,2);  y0_idx = ifind(yi,y0);
        z0 = data_in(i,l+2);
        %data_temp(y0_idx,x0_idx)=z0;
        if cust_bound
            if (x0 >= xbound(1) & x0 <= xbound(2) & y0 >= ybound(1) & y0 <= ybound(2))
                data_temp(x0_idx,y0_idx)=z0;    %non-interpolated data matrix with 99999 in masked areas
            end
        else
            data_temp(x0_idx,y0_idx)=z0;
        end
    end
      
    %now interploate
    if option
        %distances
        x=near_neighbor(xi,step_interp);
        y=near_neighbor(yi,step_interp);
        
        %data interpolated
        if sw==1
            data_temp2 = interp2(data_temp,step_interp,'nearest');
            sw=0;
        else
            data_temp2 = interp2(data_temp,step_interp,imethod);
        end

        %clear distances
        %x=xf;y=yf;
        %clear xf;clear xi;%clear yf;clear yi;

        %assign data to final structure
        data_out(:,:,l)=data_temp2;

    else %no interpolation
        data_out(:,:,l)=data_temp;
        %x=x_old;y=y_old;
        x_temp=x_old(ismember(x_old,xi));    %keeps only the elements of the bounded image, won;t affect if no bounds set
        y_temp=y_old(ismember(x_old,xi));
        y=y_temp(ismember(y_temp,yi));
        x=x_temp(ismember(y_temp,yi));
    end
    %clear data_temp;
end

%make data mask
data_mask = ones(nx0,ny0);  %make mask same size as interpolated matrix
[I,J]=find(data_out(:,:,1)==99999);
for indx=1:length(I)
   data_mask(I(indx),J(indx))=0;
end

%output  new dimensions
ux=unique(x);uy=unique(y);
if length(ux)>1
dx=abs(ux(2)-ux(1));
else
    dx=10;
end

dy=abs(uy(2)-uy(1));

%disp(sprintf('----------------------------------------------------------------------------------------------------------'))
disp(sprintf('Interp.   matrix dimensions: %3.1fmm\t (step=%2.2fmm)\t  %3.1fmm\t (step=%2.2fmm)\n',max(x)-min(x), dx,max(y)-min(y), dy))

 clear data_temp2;

