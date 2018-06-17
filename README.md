# A model of the entrance pupil and retinal landmarks
<img src="img/renderEyePose.png" height="150">

These routines implement a ray-traced model eye in MATLAB. A primary application of the model is to describe the entrance pupil in the image plane for a rotated eye. The entrance pupil is described by the parameters of an ellipse fit to the pupil perimeter, and those parameters are given in "transparent" form (center x, center y, area, non-linear eccentricity, tilt).

The model is described in:

 * GK Aguirre (2018) [The Entrance Pupil of the Human Eye](https://www.biorxiv.org/content/early/2018/05/18/325548). bioRxiv.

These routines are used to support model-based eye tracking with [transparentTrack](https://github.com/gkaguirrelab/transparentTrack)

The function `pupilProjection_fwd` implements the forward model of this projection. Inputs to this routine are:
 * `eyePose` which is a vector that describes dynamic aspects of the eye, specifically rotation in degrees of azimuth, elevation, and torsion, and the radius of the pupil aperture in mm.
 * `sceneGeometry` which is a structure that describes static aspects of the scene, including the parameters of the model eye and the properties and position of a pinhole camera. The sceneGeometry structure is generated by the function `createSceneGeometry`.

Most of the functions take key-value pairs that adjust the default behavior of the model. For example, you can generate a -3 diopter, myopic left eye that is wearing glasses with an appropriate corrective lens with the call:
```
createSceneGeometry('eyeLaterality','left','sphericalAmetropia',-3,'spectacleLens',-3)`.
```

The function `pupilProjection_inv` implements a search over eyePose parameters and executions of the forward model to find the eyePose values that best describe an observed entrance pupil ellipse.

The forward model of the appearance of the pupil and iris accounts for the refractive properties of the cornea (and any artificial lenses between the eye and the camera). The routine `virtualImageFunc.m` calculates the effect of refraction, making use of calls to `rayTraceEllipsoids.m`. An improvement in the execution time of the forward model can be achieved by compiling the ray tracing routines. To do so, issue the command `compileVirtualImageFunc` at the MATLAB console. A compiled MEX file version of `virtualImageFunc` will be placed in the `bin` directory of this repository if it is not already present. By default, the model assumes that the eye is being observed in the near infra-red range. To use refractive index values appropriate for the visible range, pass the key-value pair `'spectralDomain','vis'` when creating the scene geometry.

To install and configure the repository, first install [toolboxToolbox (tBtB)](https://github.com/ToolboxHub/ToolboxToolbox), which provides for declarative dependency management for Matlab. Once tBtB is installed, the code (and all its dependencies) will be installed and readied for use with the command `tbUse('gkaModelEye');`. If you do not wish to use tBtB, add the [quadfit toolbox](https://www.mathworks.com/matlabcentral/fileexchange/45356-fitting-quadratic-curves-and-surfaces) to your path.

A good place start is to render the model eye for different poses and examine the parameters of the pupil ellipse. This example renders an emmetropic right eye, observed in the near infra-red range, that is rotated to -30 degrees azimuth, -5 degrees elevation, and has a pupil aperture 2 mm in radius.
```
    sceneGeometry=createSceneGeometry();
    eyePose = [-30 -5 0 2];
    renderEyePose(eyePose, sceneGeometry);
    pupilEllipse = pupilProjection_fwd(eyePose,sceneGeometry);
```

The components of the model eye can be displayed in a cross-section schematic:
```
    sceneGeometry=createSceneGeometry();
    plotModelEyeSchematic(sceneGeometry.eye);
```

<p align="center">
	<img src="img/plotModelEyeSchematic.png" height="400">
</p>


A hierarchy of the functions is as follows:
```
    pupilProjection_inv
            |
            V
    pupilProjection_fwd  <--  createSceneGeometry
            |                   |-- modelEyeParameters
            V                   |     '--returnRefractiveIndex
    virtualImageFunc            |
            |                   |-- addContactLens
            V                   '-- addSpectacleLens
    rayTraceEllipsoids
```

Most functions have associated examples in the header comments. To automatically run all examples, ensure that the [ExampleTest toolbox](https://github.com/isetbio/ExampleTestToolbox.git) is on the path. This command will then test all examples:
```
	[names,status] = RunExamples(fullfile(userpath(),'toolboxes','gkaModelEye'))
```

