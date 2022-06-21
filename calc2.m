close all
model = createpde('structural','static-solid');
gm = importGeometry(model,'arch2.stl');
figure
pdegplot(model,'facelabels','on')
view(30,30);
title('Bracket with Face Labels')
% F2, F4 is fixed
structuralProperties(model,'YoungsModulus',3.25e10, ...
                           'PoissonsRatio',0.17);
structuralBC(model,'Face',1,'Constraint','fixed');
structuralBC(model,'Face',6,'Constraint','fixed');
% structuralBoundaryLoad (model,'Vertex',7,'force',[0;0;-360e3]);
% structuralBoundaryLoad (model,'Vertex',8,'force',[0;0;-360e3]);

structuralBoundaryLoad (model,'face',2,'surfacetraction',[0;0;-10.5e3*2-12e3]);
structuralBoundaryLoad (model,'face',5,'surfacetraction',[0;0;-10.5e3*2-12e3]);


mesh = generateMesh(model);
figure
pdeplot3D(model)
title('Mesh with Quadratic Tetrahedral Elements');
result = solve(model);
xx = result.Stress.xx;
yy = result.Stress.yy;
zz = result.Stress.zz;
xy = result.Stress.xy;
yz = result.Stress.yz;
xz = result.Stress.xz;

% xx = result.Stress.xx;
% xx(xx>1e5) = 1e5;
% xx(xx<-1e5) = -1e5;
% 
% yy = result.Stress.yy;
% yy(yy>1e5) = 1e5;
% yy(yy<-1e5) = -1e5;
% zz = result.Stress.zz;
% zz(zz>1e6) = 1e6;
% zz(zz<-1e6) = -1e6;

% xy = result.Stress.xy;
% xy(xy>1e5) = 1e5;
% xy(xy<-1e5) = -1e5;
% 
% yz = result.Stress.yz;
% yz(yz>1e5) = 1e5;
% yz(yz<-1e5) = -1e5;
% 
% xz = result.Stress.xz;
% xz(xz>1e5) = 1e5;
% xz(xz<-1e5) = -1e5;

figure
pdeplot3D(model,'ColorMapData',xx, 'Deformation', result.Displacement, 'mesh', 'on')
title('xx-Stress')
colormap('jet')

figure
pdeplot3D(model,'ColorMapData',(zz))
title('zz-Stress')
colormap('jet')


figure
pdeplot3D(model,'ColorMapData',(yz))
title('yz-Stress')
colormap('jet')

figure
pdeplot3D(model,'ColorMapData',(xy))
title('xy-Stress')
colormap('jet')


figure
pdeplot3D(model,'ColorMapData',(xz))
title('xz-Stress')
colormap('jet')

O = result.Stress.xx+result.Stress.yy;
R = 0.5*sqrt((result.Stress.xx-result.Stress.yy).^2 + (result.Stress.xz - result.Stress.yz).^2);
maxSigma = O+R;
minSigma = O-R;
% figure
% histogram(maxSigma)
% figure
% histogram(minSigma)


