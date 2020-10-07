function varargout = m_mollweide(varargin)
% M_MOLLWEIDE plot of data on sphere
%
% h = m_mollweide(xyzw,'Param1',Value1,'Param2',Value2,...);
% h = m_mollweide(xyz,w,'Param1',Vaule1,'Param2',Value2,...);
%
% xyz   Cartesian coordinates of vertices on sphere as either cell(nVert,1)
%       of [3,1] or [nVert,3]
% w     [nVerts,1] function values at vertices
% xyzw  Cartesian coordinates and function values as either cell(nVert,1)
%       of [x;y;z;w] or [nVert,4]
%
% h     graphics handle
%
% Parameters       Values
% 'box'            ({'on'} | 'fancy' | 'off')
% 'xTick'          number of tick marks (default 6, corresp. 60 deg)
% 'yTick'          number of tick marks (default 6, corresp. 30 deg)
% 'xTickLabels'    [label1; label2; ...]
% 'yTickLabels'    [label1; label2; ...]
% 'tickLen'        tick mark length (default 0.01)
% 'tickDir'        ({'in'} | 'out')
% 'xLabelDir'      ('end' | {'middle'})
% 'yLabelDir'      ({'end'} | 'middle')
% 'fontName'       plot title, labels, annotations font (default
%                  'Helvetica')
% 'fontSize'       default 10 (pts)
% 'fontWeight'     ({'normal'} | 'bold' | 'light' | 'demi')
% 'linewidth'      default 0.5 (pts)
% 'linestyle'      matlab line spec or 'none' (default ':' for dotted)
% 'color'          matlab color spec for lines (default 'k' for black)
% 'xLabelColor'    matlab color spec for x labels (default 'k')
% 'yLabelColor'    matlab color spec for y labels (default 'k')
% 'xAxisLocation'  ({'bottom'} | 'middle' | 'top')
% 'yAxisLocation'  ({'left'} | 'middle' | 'right')
% 'pContour'       numeric vector probability level(s) for contours
%                  (default []);
% 'lnP'            w data is log probability. Only relevant when pContour
%                  is set (true | {false})
%
% TODO: Example
% TODO: See also
% Author: Lee Samuel Finn
% Copyright 2010, 2012
% $Id: m_mollweide.m 6955 2013-10-02 12:22:26Z lsfinn $

%% Parse and validate arguments
%nargoutchk(0,1);

ip = inputParser;
ip.KeepUnmatched = true;
ip.CaseSensitive = false;

ip.addRequired('xyzw',@qXYZ);
ip.addOptional('w',[],@(w)(isnumeric(w) && all(isreal(w(:)))));
ip.addParamValue('title','',@ischar);
ip.addParamValue('pContour',[], @isnumeric);

ip.addParamValue('box','on',@(s)(any(strcmpi(s,{'on','fancy','off'}))));
ip.addParamValue('xTick',6,@isnumeric);
ip.addParamValue('yTick',6,@isnumeric);
ip.addParamValue('xTickLabels',NaN,@(c)(ischar(c)||iscell(c)));
ip.addParamValue('yTickLabels',NaN,@ischar);
ip.addParamValue('tickLen',0.01,@isnumeric);
ip.addParamValue('tickDir','in',@isnumeric);
ip.addParamValue('xLabelDir','middle',...
  @(s)(any(strcmpi(s,{'end','middle'}))));
ip.addParamValue('yLabelDir','end',@(s)(any(strcmpi(s,{'end','middle'}))));
ip.addParamValue('fontName','Helvetica',@ischar);
ip.addParamValue('fontSize',10,@isnumeric);
ip.addParamValue('fontWeight','normal',@ischar);
ip.addParamValue('linewidth',0.5,@isnumeric);
ip.addParamValue('linestyle',':',@ischar);
ip.addParamValue('color','k',@ischar);
ip.addParamValue('xLabelColor','k',@ischar);
ip.addParamValue('yLabelColor','k',@ischar);
ip.addParamValue('xAxisLocation','middle',...
  @(s)(any(strcmpi(s,{'bottom','middle','top'}))));
ip.addParamValue('yAxisLocation','left',...
  @(s)(any(strcmpi(s,{'left','middle','right'}))));
ip.addParamValue('lnP',false,@islogical);

ip.parse(varargin{:});

%%% Convert input to canonical form
% xyz = [nVert,3];
% w = [nVert,1];
tmp = ip.Results.xyzw';
if isnumeric(tmp)
  xyz = tmp(1:3,:)';
else
  xyz = cellfun(@(x)([x(1),x(2),x(3)]),tmp,'UniformOutput',false);
  xyz = cell2mat(reshape(xyz,[],1));
end

w = ip.Results.w;
if any(strcmpi('w',ip.UsingDefaults))
  if isnumeric(tmp)
    w = tmp(4,:)';
  else
    w = reshape(cellfun(@(w)(w(4)),tmp),[],1);
  end
elseif isnumeric(w)
  w = reshape(w,[],1);
else
  w = reshape(cell2mat(w),[],1);
end

%%% Initialize constants
msgIdDelDup = 'MATLAB:DelaunayTri:DupPtsWarnId';
msgIdTriDup = 'MATLAB:TriScatteredInterp:DupPtsAvValuesWarnId';
wsDel = warning('QUERY',msgIdDelDup);
wsTri = warning('QUERY',msgIdTriDup);
warning('OFF',msgIdDelDup);
warning('OFF',msgIdTriDup);

%% Build patches
% The Delaunay triangulation of the user provided vertex data has a convex
% hull whose simplices we take to be our patches. The projection of
% spherical data involves either a puncture of the sphere, or a cut through
% the sphere. Patches that cross the cut, or include the punture, must be
% treated specially.
%
% For projections involving a puncture we augment the list of vertices with
% the location of the puncture, using bilinar interpolation of the
% enclosing patch to find the function value at the new vertex.
%
% For projections involving a cut we augment the list of vertices with new
% vertices located at every location where the great circle defining a
% patch edge crosses the cut. The function value at the new vertex is
% obtained by linear interpolation along the great circle.
%
% * Add poles to list of vertices
% * Find convex hull of Delaunay triangulation of vertices
% * Find and add new vertices
% * Re-evaluate convex hull
% * Form patches from convex hull simplices
% * Plot the patches

%%% Form Delaunay Triangulation of vertices
dt = DelaunayTri(xyz);

%%% Add polar vertices
zmax = max(dt.X(:,3));
% add polar vertex if not already present
if 1-zmax > eps(1)
  % Find simplex including pole
  % Find vertices of triangular simplex face that includes pole
  % Find pole as linear combination of simplex vertices
  % Find value at pole as linear combination of simplex vertex values
  % Add new vertex to triangulation
  [si,bc] = pointLocation(dt,[0,0,zmax]);
  np = dot(w(dt.Triangulation(si,:)),bc);
  % When adding new vertex make sure it is not rejected as a duplicate
  lastwarn('');
  dt.X(end+1,:) = [0,0,1];
  [~,msgid] = lastwarn;
  if ~any(strcmpi(msgid,{msgIdDelDup,msgIdTriDup}))
    w(end+1) = np;
  end
end

zmin = min(dt.X(:,3));
if abs(1+zmin) > eps(1)
  si = pointLocation(dt,[0,0,zmin]);
  np = dot(w(dt.Triangulation(si,:)),bc);
  % When adding new vertex make sure it is not rejected as a duplicate
  lastwarn('');
  dt.X(end+1,:) = [0,0,-1];
  [~,msgid] = lastwarn;
  if ~any(strcmpi(msgid,{msgIdDelDup,msgIdTriDup}))
    w(end+1) = np;
  end
end

%%% Form triangulation's convex hull
% ch is [nSmp,nDim] array of vertex indices, with
%   nSmp is the number of simplices in the convex hull
%   nDim is the number of vertices in the each simplex
% The vertices in simplex k are indexed by ch(k,:)
ch = convexHull(dt);
nSmp = size(ch,1);

%%% Identify new vertices and add to triangulation
for kSmp = 1:nSmp
  % Get simplex vertices and values at vertices
  ndx = ch(kSmp,:);
  abc = dt.X(ndx,:);
  f = w(ndx);
  
  % Find simplex edges
  edges = diff([abc;abc(1,:)],1,1);
  da = edges(:,1);
  db = edges(:,2);
  dc = edges(:,3);
  
  % Find distance along edge to ac plane
  r = -abc(:,2)./db;
  a = abc(:,1) + r.*da;
  
  % r = 0, 1 are the coordinates along the edge of the two vertices. If the
  % intercept is in (0,1) and a < 0 then the simplex edge crosses the
  % dateline
  mask = 0 < r & r < 1 & a < 0;
  
  if any(mask)
    % Simplex edge crosses dateline
    
    % Find z-coordinate of crossing
    c = abc(:,3) + r.*dc;
    
    % New vertices
    nVert = zeros(sum(mask),3);
    nVert(:,1) = reshape(a(mask),[],1);
    nVert(:,3) = reshape(c(mask),[],1);
    for k = 1:sum(mask)
      nVert(k,:) = nVert(k,:)/norm(nVert(k,:));
    end
    
    % Interpolate function to new vertices
    abc4 = [abc;abc(1,:)];
    f4 = [f;f(1)];
    qABC4 = quaternion(zeros(4,1),abc4(:,1),abc4(:,2),abc4(:,3));
    qNVert = quaternion(...
      zeros(sum(mask),1),...
      nVert(:,1),nVert(:,2),nVert(:,3));
    kVert = 0;
    fN = zeros(1,sum(mask));
    for k = 1:3
      if mask(k)
        kVert = kVert+1;
        qA = qABC4(k);
        qB = qABC4(k+1);
        qV = qNVert(kVert);
        if abs(qV-qA) > 2*eps
          r = abs(log(qV/qA)/log(qB/qA));
        else
          r = 0;
        end
        fN(kVert) = f4(k)+(f4(k+1)-f4(k))*r;
      end
    end
    
    for k = 1:sum(mask)
      % When adding new vertex make sure it is not rejected as a duplicate
      lastwarn('');
      dt.X(end+1,:) = nVert(k,:);
      [~,msgid] = lastwarn;
      if ~any(strcmpi(msgid,{msgIdDelDup,msgIdTriDup}))
        w(end+1) = fN(k);  %#ok
      end
    end
  end
end

%% Form and plot patches
% Re-evaluate convex hull
ch = convexHull(dt);

% Convex hull simplices
x = reshape(dt.X(ch,1),size(ch));
y = reshape(dt.X(ch,2),size(ch));
z = reshape(dt.X(ch,3),size(ch));

% Naive conversion of simplex vertices to longitude, latitude
ln = atan2(y,x)*180/pi;
lt = asin(z)*180/pi;

% Address dateline issues
% all simplices are triangles
% sign of 180 should be sign of not 180
for ks = 1:length(ln)
  m180 = abs(abs(ln(ks,:))-180)<2*eps;
  if any(m180)
    % at least one vertex on dateline
    % get sign(longitude) of other vertices
    sn180 = sign(sum(sign(ln(ks,~m180))));
    % set sign(longitude) for dateline vertices appropriately
    ln(ks,m180) = 180*sn180;
  end
end

f = reshape(w(ch),size(ch))';
ln = ln';
lt = lt';

% Start figure
m_proj('mollweide','lat',[-90 90],'lon',[-180,180]);
gridOpts = rmfield(ip.Results,...
  {'xyzw','w','pContour','title','lnP'});
% if ~isfield(gridOpts,'getSet')
%   gridOpts.getSet = 'none';
% end
gridArgs = struct2cell(gridOpts);
gridArgsOpts = cell(2,numel(gridArgs));
gridArgsOpts(1,:) = fieldnames(gridOpts);
gridArgsOpts(2,:) = gridArgs;
m_grid(gridArgsOpts{:});
[xx,yy] = m_ll2xy(ln,lt);
%h = patch(xx,yy,f,'EdgeColor','none','FaceColor','interp');
h = patch(xx,yy,f,'EdgeColor','none','FaceColor','flat');
% axis equal
% axis off;

%%% Probability contours
if ~isempty(ip.Results.pContour)
  % Assume w are probabilities and that vertices are equispaced in area on
  % sphere. Plot smallest area contours enclosing pContour fraction of
  % probability.
  
  %%% Create interpolating function
  phix = [atan2(dt.X(:,2),dt.X(:,1)), asin(dt.X(:,3))];
  
  % extend domain for purposes of interpolation
  phix = [...
    bsxfun(@plus,phix,[-2*pi,0]); ...
    phix; ...
    bsxfun(@plus,phix,[2*pi,0])];
  f = TriScatteredInterp(phix,[w;w;w],'natural');
  
  %%% Interpolate onto mesh
  nMesh = 2*ceil(sqrt(numel(w)));
  x = (-nMesh:nMesh)/(nMesh+1);
  phi = x*pi;
  [phi,x] = meshgrid(phi,x);
  mw = f(phi,x);
  
  % Sort pixels from lowest to highest probability and find cumulative
  % probability
  pk = sort(mw(:),'ascend');
  % if lnP true then w is log(p)
  if ip.Results.lnP
    cpk = cumsum(exp(pk));
  else
    cpk = cumsum(pk);
  end
  cpk = cpk/cpk(end);
  
  % Find the probability level
  prob = ip.Results.pContour;
  cLvl = zeros(numel(prob),1);
  for kLvl = 1:numel(prob)
    pix = find(cpk<(1-prob(kLvl)),1,'last');
    r = ((1-prob(kLvl))-cpk(pix))/(cpk(pix+1)-cpk(pix));
    cLvl(kLvl) = pk(pix) + r*(pk(pix+1)-pk(pix));
    nPix = numel(pk)-pix;
    aPix = (4*pi/numel(pk))*(180/pi)^2;
    fprintf('%5.1f%% Contour = %-7.3g deg^2 (nPix = %d)\n',...
      100*prob(kLvl),nPix*aPix,nPix);
  end
  
  %%% draw contour
  hold on;
  m_contour(phi*180/pi,asin(x)*180/pi,mw,cLvl,'-k','linewidth',2);
  hold off;
  
  %   % For contouring a single level tricontour requires that the level be
  %   % listed twice
  %   if 1 == numel(prob)
  %     cLvl = cLvl*[1,1];
  %   end
  %
  %   % Contour directly on projection
  %   lt = asin(dt.X(:,3))*180/pi;
  %   ln = atan2(dt.X(:,2),dt.X(:,1))*180/pi;
  %   [xx,yy] = m_ll2xy(ln,lt);
  %   hold on;
  %   tricontour([xx(:),yy(:)],ch,w,cLvl,'EdgeColor','k');
  %   hold off;
  
end

%% Colorbar, lighting
if isfield(gridOpts,'colorbar') && ~ischar(gridOpts.colorbar)
  colorbar(gridOpts.colorbar);
end
colormap(jet(512));
lighting phong;

%% Add title
if ~isempty(ip.Results.title)
  title(ip.Results.title);
end


if 1 == nargout
  varargout{1} = h;
end

warning(wsDel);
warning(wsTri);

return


function tf = qXYZ(x)
% QXYZ validate vertex list input argument
% Returns true if x is either
%   cell(nVert,1) of [3,1] or [4,1]
% or
%   [nVert,3] or [nVert,4]

if iscell(x) && length(x) == numel(x)
  n = cellfun('prodofsize',x);
  tf = all(4 == n) || all(3 == n);
elseif isnumeric(x)
  n = size(x,2);
  tf = 3 == n || 4 == n;
else
  tf = false;
end
return
