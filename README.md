# Forms: Flexible Object Reconstruction from Multiple Silhouettes

### Thomas J. Cashman, University of Lugano
### Andrew W. Fitzgibbon, Microsoft Research, Cambridge

### Updates by Jason M. Manley, Rockefeller University
### Questions? jmanley@rockefeller.edu

See the bottom for the original documentation provided by Cashman & Fitzgibbon
to reproduce the results from their 2012 paper (http://ieeexplore.ieee.org/document/6165306/).
The dolphin model produced by running the code described in their README.txt is given in _fitzgibbonModel.mat_.

This code was written and tested in *MATLAB R2017a* and requires the *Optimization Toolbox*.

Original source code retrieved from https://archive.codeplex.com/?p=forms.


## What do I do?

*Forms* uses a procedure developed by Cashman & Fitzgibbon to fit a 3D morphable
model to a collection of 2D images showing some type of object in its various poses
and conformations. This fitting is formulated as an energy minimization problem,
and the model represents the object as a linear combination of subdivision surfaces.


## Summary of updates from the original source code

- _set_initial_conditions.m_ GUI for easily inputting constraints and initial orientations for new datasets.
- Reformatting the _forms.m_ output as a structure (instead of a vector with very careful indexing).
- Once a model has been built, _fit_model_to_image.m_ allows fitting a model to another image outside the original input set.
- _calculate_shapediff.m_ provides a potential metric for how well the model has fit the image's silhouette.


## Example scripts

See _run_forms_bananas_dolphins.m_ for fitting the model to the example data provided by Cashman & Fitzgibbon.

See _run_forms_new_images.m_ for curating new images and fitting the model to a new dataset.


## How to structure input data

This algorithm requires a set of 2D images, in which the desired object has
been segmented. While Cashman & Fitzgibbon stored their data in .ply and .fpj
files (as seen in forms/projects), we have opted to simply store our data in
structure arrays within .mat files. This is referred to as a Forms 'project'.


#### Each project should contain the following structure fields:

_project.vertices_    -> a Nx3 matrix specifying the N vertices of the pre-defined
                         mesh model

_project.mesh_        -> a meshtri (see Cashman & Fitzgibbon's class definition)

_project.faces_       -> a Nfx3 matrix specifying the indices of vertices which
                         form the Nf triangular faces of the pre-defined mesh model

_project.cand_ixs_    ->

_project.cand_uvs_    ->

_project.cand_limits_ ->

_project.cand_dists_  ->

_project.cand_derivs_ ->


For our analyses to date, the above fields have been the same as those utilized
in Cashman & Fitzgibbon's dolphin analysis (see _fitzgibbonModel.mat_). The following
fields should be specified for future analyses of new dolphins datasets:

_project.images_      -> a Nx1 structure array containing the input dataset for
                         N input images

_project.parameters_  -> a structure array containing the parameters for running
                         the model fitting in forms.m. See makeModelParameters.m.


_project.images_ should contain the following fields for the i-th image project.images(i):


#### The following fields must be directly specified by the user.

_image_            -> the mxn(x3) original image

_silhouette_       -> a Nx2 matrix containing the 2D location of N points along
                      the object's silhouette

_frame_            -> the frame index at which this image occurs in the movie

_movie_            -> the path to this frame's original movie


#### The following fields will be found using _set_initial_conditions.m_.

_constraints3d_    -> a Nx1 vector containing the indices of the N contraint points
                      on the 3D mesh model. This index refers to an actual 3D location
                      given in project.vertices.

_constraints2d_    -> a Nx2 matrix containing the 2D locations of the N constraint
                      points on the image, such that constraints3d(i) corresponds to
                      the location of constraints2d(i,:) on the model.

_constraintsonsil_ -> a Nx1 binary vector describing whether is constraint point is
                      on (1) or off (0) the silhouette.

_normalsLeft_      -> a Nx1 binary vector describing whether the normal vectors
                      should be flipped (0) or not (1).

_rotate_           -> a rotation matrix describing the approximate orientation
                      of the dolphin for the model's initial conditions. This is
                      specified by hand using set_initial_conditions.m.


#### The following fields will be automatically found using the pipeline below.

_translate_        -> a translation matrix for aligning the 3D model with the
                      2D silhouette, found using align_project.m.

_scale_            -> a dilation matrix for aligning the 3D model with the
                      2D silhouette, found using align_project.m.

_transform_        -> the full transformation matrix, such that
                      transform = rotate * translate * scale.

_points_           -> a 4x(N/4)x2 matrix containing a reshaped version of the
                      silhouette, in order to be compatible with forms.m. % Why?


## How to fit a model to the input data

Once you have initialized a project with a mesh structure, images, and silhouettes,
follow the pipeline below to find a 3D morphable model for your object.

Firstly, the initial conditions must be specified by running the following
for each of the images i=1:length(project.images):

> set_initial_conditions(project,i)

This will open a figure and guide you through setting the initial orientation
of the object (using the angle sliders), selecting the constraint points, and
properly orienting the silhouette normal vectors. Note: this also automatically
runs the alignment code > align_project(project).

If you choose to use the default parameters, the following creates the
parameter structure:

> parameters = makeModelParameters;   % see makeModelParameters or Cashman & Fitzgibbon's original paper for a description of these parameters

The model can then by found by running:

> model = forms(project, parameters);

The fit of the model for image i can the be visualized by running:

> plot_modelfit(project, i, model);

The norm of the ith basis shape can be plotted on the 0th mode by running:

> plot_basisshapenorm(project, model, i);


See Cashman & Fitzgibbon's original README below. Note that the forms
function has been modified in order to store the model in a structure, rather
than in a carefully indexed single vector. The function _archive/convert_modelVec_to_modelStruct_
can be utilized to convert the vector version of a model to the structure
version of a model that is utilized here. In addition, the parameters inputs
have been modified to be stored in a structure.


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Forms: Flexible Object Reconstruction from Multiple Silhouettes

Thomas J. Cashman          Andrew W. Fitzgibbon
University of Lugano       Microsoft Research, Cambridge


This MATLAB code gives the implementation for our paper 'What Shape are
Dolphins? Building 3D Morphable Models from 2D Images'. You will need the
MATLAB Optimization Toolbox to run it.


To try the code on one of our datasets, load a Forms 'project', and then use the
function 'forms' to run our optimization. For example:

>> bananas = read_project('.\projects\bananas.fpj');
>> banana_model = forms(bananas, 2);

You can then compare the fit of the model before and after optimization, by
using

>> plot_modelfit(bananas, 1);                % Show banana 1 before optimization
>> plot_modelfit(bananas, 1, banana_model);  % Show banana 1 after optimization

Or take a look at the basis shapes plotted as a colour map by using

>> plot_basisshapenorm(bananas, banana_model, 1);  % First basis shape norm
>> plot_basisshapenorm(bananas, banana_model, 2);  % Second basis shape norm


To reproduce our results, use the weights and normal noise estimates described
in the paper, and use the function 'align_project' to find the camera
translation and scale parameters automatically, i.e.

>> bananas = read_project('.\projects\bananas.fpj');
>> bananas = align_project(bananas);
>> banana_model = forms(bananas, 2, 0.5, 0.25, 4);

>> pigeons = read_project('.\projects\pigeons.fpj');
>> pigeons = align_project(pigeons);
>> pigeon_model = forms(pigeons, 7, 0.25, 0.05, 5);

>> bears = read_project('.\projects\bears.fpj');
>> bears = align_project(bears);
>> bear_model = forms(bears, 10, 0.25, 0.25, 4);

Note that you should expect some of these optimizations to take a long time.
This is particularly true for the dolphins project, which can be used to build a
a model from 32 dolphin instances:

>> dolphins = read_project('.\projects\dolphins.fpj');
>> dolphins = align_project(dolphins);
>> dolphin_model = forms(dolphins, 8, 0.5, 0.25, 10);

The dolphin dataset is available here with manual segmentations, rather than the
grab cut silhouettes from Powerpoint 2010. This is the reason we suggest
sigma_norm = 10 above, rather than the value 40 / 3 that appears in the paper:
the manual segmentations have more reliable normals, so we can use this slightly
lower noise estimate. However, there is no significant difference between the
manual and automatic segmentations; the results are very similar.
