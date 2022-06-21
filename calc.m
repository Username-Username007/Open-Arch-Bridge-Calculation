close all
model = createpde('structural','static-solid');
gm = importGeometry(model,'arch.stl');
figure
pdegplot(model,'FaceLabels','on')
view(30,30);
title('Bracket with Face Labels')
% F2, F4 is fixed
structuralProperties(model,'YoungsModulus',3.25e10, ...
                           'PoissonsRatio',0.17);
structuralBC(model,'Face',4,'Constraint','fixed');
structuralBC(model,'Face',2,'Constraint','fixed');
structuralBoundaryLoad (model,'Face',3,'SurfaceTraction',[0;-10.5e3;0]);
mesh = generateMesh(model);
figure
pdeplot3D(model)
title('Mesh with Quadratic Tetrahedral Elements');
result = solve(model);

figure
pdeplot3D(model,'ColorMapData',result.Stress.xx, 'Deformation', result.Displacement, 'mesh', 'on')
title('xx-Stress')
colormap('jet')
x
figure
pdeplot3D(model,'ColorMapData',result.Stress.yy)
title('yy-Stress')
colormap('jet')

O = result.Stress.xx+result.Stress.yy;
R = 0.5*sqrt((result.Stress.xx-result.Stress.yy).^2 + (result.Stress.xz - result.Stress.yz).^2);
maxSigma = O+R;
minSigma = O-R;
figure
histogram(maxSigma)
figure
histogram(minSigma)


