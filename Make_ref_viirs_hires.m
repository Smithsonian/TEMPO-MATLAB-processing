function VIIRS = Make_ref_viirs_hires(fn_dnb,LAT0,LON0,up_sampling)

% Reference VIIRS-DNB Composite
dnb = imread(fn_dnb);
dnb(dnb<1) = 1;
dnb(dnb>1e4) = 1e4;

% Upsampled site indices
ii = ((0.5+1/(2*up_sampling(2))):(1/up_sampling(2)):(size(LAT0,1)+0.5-1/(2*up_sampling(2))))';
jj = (0.5+1/(2*up_sampling(1))):(1/up_sampling(1)):(size(LAT0,2)+0.5-1/(2*up_sampling(1)));
II = repmat(ii,1,length(jj));
JJ = repmat(jj,length(ii),1);

LAT = interp2(LAT0,JJ,II);
LON = interp2(LON0,JJ,II);

Y = (75-LAT)*3600/15+1;
X = (LON+180)*3600/15+1;

% figure, plot(X(:),Y(:),'.')

X(X<1) = 1;
Y(Y<1) = 1;
X(X>size(dnb,2)) = size(dnb,2);
Y(Y>size(dnb,1)) = size(dnb,1);

VIIRS = interp2(1:size(dnb,2),1:size(dnb,1),dnb,X,Y);

% figure, imagesc(fliplr(CITY),[0 100]), colormap turbo
% figure, imagesc(fliplr(VIIRS),[0 100]), colormap turbo

% compo = zeros(size(CITY,1),size(CITY,2),3);
% compo(:,:,1) = fliplr(CITY);
% compo(:,:,2) = fliplr(VIIRS);
% compo(:,:,3) = fliplr(CITY);
% 
% compo(compo<0) = 0;
% compo(compo>100) = 100;