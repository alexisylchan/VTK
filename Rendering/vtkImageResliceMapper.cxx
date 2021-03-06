/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageResliceMapper.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImageResliceMapper.h"

#include "vtkRenderer.h"
#include "vtkCamera.h"
#include "vtkImageSlice.h"
#include "vtkImageData.h"
#include "vtkImageProperty.h"
#include "vtkLookupTable.h"
#include "vtkMath.h"
#include "vtkMatrix4x4.h"
#include "vtkAbstractTransform.h"
#include "vtkPlane.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkGarbageCollector.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkImageResliceToColors.h"
#include "vtkGraphicsFactory.h"

//----------------------------------------------------------------------------
// Needed when we don't use the vtkStandardNewMacro.
vtkInstantiatorNewMacro(vtkImageResliceMapper);

//----------------------------------------------------------------------------
vtkImageResliceMapper* vtkImageResliceMapper::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkGraphicsFactory::CreateInstance("vtkImageResliceMapper");
  return static_cast<vtkImageResliceMapper *>(ret);
}

//----------------------------------------------------------------------------
vtkImageResliceMapper::vtkImageResliceMapper()
{
  this->ImageReslice = vtkImageResliceToColors::New();
  this->ResliceMatrix = vtkMatrix4x4::New();
  this->WorldToDataMatrix = vtkMatrix4x4::New();
  this->SliceToWorldMatrix = vtkMatrix4x4::New();
}

//----------------------------------------------------------------------------
vtkImageResliceMapper::~vtkImageResliceMapper()
{
  if (this->ImageReslice)
    {
    this->ImageReslice->Delete();
    }
  if (this->ResliceMatrix)
    {
    this->ResliceMatrix->Delete();
    }
  if (this->WorldToDataMatrix)
    {
    this->WorldToDataMatrix->Delete();
    }
  if (this->SliceToWorldMatrix)
    {
    this->SliceToWorldMatrix->Delete();
    }
}

//----------------------------------------------------------------------------
void vtkImageResliceMapper::SetSlicePlane(vtkPlane *plane)
{
  if (this->SlicePlane == plane)
    {
    return;
    }
  if (this->SlicePlane)
    {
    this->SlicePlane->Delete();
    }
  if (!plane)
    {
    this->SlicePlane = vtkPlane::New();
    }
  else
    {
    this->SlicePlane = plane;
    plane->Register(this);
    }

  this->Modified();
}

//----------------------------------------------------------------------------
void vtkImageResliceMapper::ReleaseGraphicsResources(vtkWindow *)
{
  // see OpenGL subclass for implementation
}

//----------------------------------------------------------------------------
void vtkImageResliceMapper::Render(vtkRenderer *, vtkImageSlice *)
{
  // see OpenGL subclass for implementation
}

//----------------------------------------------------------------------------
int vtkImageResliceMapper::ProcessRequest(
  vtkInformation* request, vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  if (request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_INFORMATION()))
    {
    // Get point/normal from camera
    if (this->SliceFacesCamera || this->SliceAtFocalPoint)
      {
      vtkRenderer *ren = this->GetCurrentRenderer();

      if (ren)
        {
        vtkCamera *camera = ren->GetActiveCamera();

        if (this->SliceAtFocalPoint)
          {
          this->SlicePlane->SetOrigin(camera->GetFocalPoint());
          }
        if (this->SliceFacesCamera)
          {
          double normal[3];
          camera->GetDirectionOfProjection(normal);
          normal[0] = -normal[0];
          normal[1] = -normal[1];
          normal[2] = -normal[2];
          this->SlicePlane->SetNormal(normal);
          }
        }
      }

    // use superclass method to update other important info
    return this->Superclass::ProcessRequest(
      request, inputVector, outputVector);
    }

  return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}

//----------------------------------------------------------------------------
// Update the WorldToData transformation matrix, which is just the
// inverse of the vtkProp3D matrix.
void vtkImageResliceMapper::UpdateWorldToDataMatrix(vtkImageSlice *prop)
{
  // copy the matrix, but only if it has changed (we do this to
  // preserve the modified time of the matrix)
  double tmpmat[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };
  if (!prop->GetIsIdentity())
    {
    vtkMatrix4x4::Invert(*prop->GetMatrix()->Element, tmpmat);
    }
  double *mat = *this->WorldToDataMatrix->Element;
  for (int i = 0; i < 16; i++)
    {
    if (mat[i] != tmpmat[i])
      {
      this->WorldToDataMatrix->DeepCopy(tmpmat);
      break;
      }
    }
}

//----------------------------------------------------------------------------
// Update the SliceToWorld transformation matrix
void vtkImageResliceMapper::UpdateSliceToWorldMatrix(vtkCamera *camera)
{
  // Get slice plane in world coords by passing null as the prop matrix 
  double plane[4];
  this->GetSlicePlaneInDataCoords(0, plane);

  // Make sure normal is facing towards camera
  vtkMatrix4x4 *viewMatrix = camera->GetViewTransformMatrix();
  double *ndop = viewMatrix->Element[2];
  if (vtkMath::Dot(ndop, plane) < 0)
    {
    plane[0] = -plane[0];
    plane[1] = -plane[1];
    plane[2] = -plane[2];
    plane[3] = -plane[3];
    }

  // The normal is the first three elements
  double *normal = plane;

  // The last element is -dot(normal, origin)
  double dp = -plane[3];

  // Compute rotation angle between camera axis and slice normal
  double vec[3];
  vtkMath::Cross(ndop, normal, vec);
  double costheta = vtkMath::Dot(ndop, normal);
  double sintheta = vtkMath::Norm(vec);
  double theta = atan2(sintheta, costheta);
  if (sintheta != 0)
    {
    vec[0] /= sintheta;
    vec[1] /= sintheta;
    vec[2] /= sintheta;
    }
  // convert to quaternion
  costheta = cos(0.5*theta);
  sintheta = sin(0.5*theta);
  double quat[4];
  quat[0] = costheta;
  quat[1] = vec[0]*sintheta;
  quat[2] = vec[1]*sintheta;
  quat[3] = vec[2]*sintheta;
  // convert to matrix
  double mat[3][3];
  vtkMath::QuaternionToMatrix3x3(quat, mat);

  // Create a slice-to-world transform matrix
  // The columns are v1, v2, normal
  vtkMatrix4x4 *sliceToWorld = this->SliceToWorldMatrix;

  double v1[3], v2[3];
  vtkMath::Multiply3x3(mat, viewMatrix->Element[0], v1);
  vtkMath::Multiply3x3(mat, viewMatrix->Element[1], v2);
  vtkMath::Multiply3x3(mat, viewMatrix->Element[2], vec);

  sliceToWorld->Element[0][0] = v1[0];
  sliceToWorld->Element[1][0] = v1[1];
  sliceToWorld->Element[2][0] = v1[2];
  sliceToWorld->Element[3][0] = 0.0;

  sliceToWorld->Element[0][1] = v2[0];
  sliceToWorld->Element[1][1] = v2[1];
  sliceToWorld->Element[2][1] = v2[2];
  sliceToWorld->Element[3][1] = 0.0;

  sliceToWorld->Element[0][2] = normal[0];
  sliceToWorld->Element[1][2] = normal[1];
  sliceToWorld->Element[2][2] = normal[2];
  sliceToWorld->Element[3][2] = 0.0;

  sliceToWorld->Element[0][3] = -dp*normal[0];
  sliceToWorld->Element[1][3] = -dp*normal[1];
  sliceToWorld->Element[2][3] = dp-dp*normal[2];
  sliceToWorld->Element[3][3] = 1.0;
}

//----------------------------------------------------------------------------
// Do all the fancy math to set up the reslicing
void vtkImageResliceMapper::UpdateResliceInformation(vtkRenderer *ren)
{
  vtkMatrix4x4 *resliceMatrix = this->ResliceMatrix;
  vtkImageResliceToColors *reslice = this->ImageReslice;

  // Create the reslice matrix by multiplying by the prop's matrix
  vtkMatrix4x4::Multiply4x4(
    this->WorldToDataMatrix, this->SliceToWorldMatrix, resliceMatrix);

  // Get the projection matrix
  double aspect = ren->GetTiledAspectRatio();
  vtkCamera *camera = ren->GetActiveCamera();
  vtkMatrix4x4 *viewMatrix = camera->GetViewTransformMatrix();
  vtkMatrix4x4 *projMatrix = camera->GetProjectionTransformMatrix(
                               aspect, 0, 1);

  // Compute other useful matrices
  double worldToView[16];
  double viewToWorld[16];
  double planeWorldToView[16];
  vtkMatrix4x4::Multiply4x4(
    *projMatrix->Element, *viewMatrix->Element, worldToView);
  vtkMatrix4x4::Invert(worldToView, viewToWorld);
  vtkMatrix4x4::Transpose(viewToWorld, planeWorldToView);

  double worldToSlice[16];
  double viewToSlice[16];
  vtkMatrix4x4::Invert(*this->SliceToWorldMatrix->Element, worldToSlice);
  vtkMatrix4x4::Multiply4x4(worldToSlice, viewToWorld, viewToSlice);

  // Get slice plane in world coords by passing null as the matrix
  double plane[4];
  this->GetSlicePlaneInDataCoords(0, plane);

  // Check whether normal is facing towards camera, the "ndop" is
  // the negative of the direction of projection for the camera
  double *ndop = viewMatrix->Element[2];
  if (vtkMath::Dot(ndop, plane) < 0)
    {
    plane[0] = -plane[0];
    plane[1] = -plane[1];
    plane[2] = -plane[2];
    plane[3] = -plane[3];
    }

  // Get the z position of the slice in slice coords
  // (requires plane to be normalized by GetSlicePlaneInDataCoords)
  double z = (plane[2] - 2.0)*plane[3];

  // Transform the plane into view coordinates, using the transpose
  // of the inverse of the world-to-view matrix
  vtkMatrix4x4::MultiplyPoint(planeWorldToView, plane, plane);

  // Compute the bounds in slice coords
  double xmin = VTK_DOUBLE_MAX;
  double xmax = -VTK_DOUBLE_MAX;
  double ymin = VTK_DOUBLE_MAX;
  double ymax = -VTK_DOUBLE_MAX;

  for (int i = 0; i < 4; i++)
    {
    // The four corners of the view
    double x = (((i & 1) == 0) ? -1.0 : 1.0);
    double y = (((i & 2) == 0) ? -1.0 : 1.0);

    double hpoint[4];
    hpoint[0] = x;
    hpoint[1] = y;
    hpoint[2] = 0.0;
    hpoint[3] = 1.0;

    if (fabs(plane[2]) < 1e-6)
      {
      // Looking at plane edge-on, just put some
      // points at front clipping plane, others at back plane
      hpoint[2] = (((i & 1) == 0) ? 0.0 : 1.0);
      }
    else
      {
      // Intersect with the slice plane
      hpoint[2] = - (x*plane[0] + y*plane[1] + plane[3])/plane[2];

      // Clip to the front and back clipping planes
      if (hpoint[2] < 0)
        {
        hpoint[2] = 0.0;
        }
      else if (hpoint[2] > 1)
        {
        hpoint[2] = 1.0;
        }
      }

    // Transform into slice coords
    vtkMatrix4x4::MultiplyPoint(viewToSlice, hpoint, hpoint);

    x = hpoint[0]/hpoint[3];
    y = hpoint[1]/hpoint[3];

    // Find min/max in slice coords
    if (x < xmin) { xmin = x; }
    if (x > xmax) { xmax = x; }
    if (y < ymin) { ymin = y; }
    if (y > ymax) { ymax = y; }
    }

  // The ResliceExtent is always set to the renderer size,
  // this is the maximum size ever required and sticking to
  // this size avoids any memory reallocation on GPU or CPU
  int *size = ren->GetSize();
  int xsize = ((size[0] <= 0) ? 1 : size[0]);
  int ysize = ((size[1] <= 0) ? 1 : size[1]);

  int extent[6];
  extent[0] = 0;
  extent[1] = xsize - 1;
  extent[2] = 0;
  extent[3] = ysize - 1;
  extent[4] = 0;
  extent[5] = 0;

  // Find the spacing
  double spacing[3];
  spacing[0] = (xmax - xmin)/xsize;
  spacing[1] = (ymax - ymin)/ysize;
  spacing[2] = 1.0;

  // Corner of resliced plane, including half-pixel offset to
  // exactly match texels to pixels in the final rendering
  double origin[3];
  origin[0] = xmin + 0.5*spacing[0];
  origin[1] = ymin + 0.5*spacing[1];
  origin[2] = z;

  // Prepare for reslicing
  reslice->SetResliceAxes(resliceMatrix);
  reslice->SetOutputExtent(extent);
  reslice->SetOutputSpacing(spacing);
  reslice->SetOutputOrigin(origin);

  if (this->SliceFacesCamera)
    {
    // if slice follows camera, use reslice to set the border
    reslice->SetBorder(this->Border);
    }
  else
    {
    // tell reslice to use a double-thickness border,
    // since the polygon geometry will dictate the actual size
    reslice->SetBorder(2);
    }
}

//----------------------------------------------------------------------------
// Do all the fancy math to set up the reslicing
void vtkImageResliceMapper::UpdateColorInformation(vtkImageProperty *property)
{
  vtkScalarsToColors *lookupTable = this->DefaultLookupTable;

  if (property)
    {
    double colorWindow = property->GetColorWindow();
    double colorLevel = property->GetColorLevel();
    if (property->GetLookupTable())
      {
      lookupTable = property->GetLookupTable();
      if (!property->GetUseLookupTableScalarRange())
        {
        lookupTable->SetRange(colorLevel - 0.5*colorWindow,
                              colorLevel + 0.5*colorWindow);
        }
      }
    else
      {
      lookupTable->SetRange(colorLevel - 0.5*colorWindow,
                            colorLevel + 0.5*colorWindow);
      }
    }
  else
    {
    lookupTable->SetRange(0, 255);
    }

  this->ImageReslice->SetLookupTable(lookupTable);
}

//----------------------------------------------------------------------------
// Do all the fancy math to set up the reslicing
void vtkImageResliceMapper::UpdateResliceInterpolation(
  vtkImageProperty *property)
{
  // set the interpolation mode and border
  int interpMode = VTK_RESLICE_NEAREST;

  if (property)
    {
    switch(property->GetInterpolationType())
      {
      case VTK_NEAREST_INTERPOLATION:
        interpMode = VTK_RESLICE_NEAREST;
        break;
      case VTK_LINEAR_INTERPOLATION:
        interpMode = VTK_RESLICE_LINEAR;
        break;
      case VTK_CUBIC_INTERPOLATION:
        interpMode = VTK_RESLICE_CUBIC;
        break;
      }
    }

  this->ImageReslice->SetInterpolationMode(interpMode);
}

//----------------------------------------------------------------------------
void vtkImageResliceMapper::CheckerboardImage(
  vtkImageData *input, vtkCamera *camera, vtkImageProperty *property)
{
  // Use focal point as center of checkerboard pattern.  This guarantees
  // exactly the same checkerboard for all images in the scene, which is
  // useful when doing multiple overlays.
  double focalPoint[4];
  camera->GetFocalPoint(focalPoint);
  focalPoint[3] = 1.0;

  double worldToSlice[16];
  vtkMatrix4x4::Invert(*this->SliceToWorldMatrix->Element, worldToSlice);

  vtkMatrix4x4::MultiplyPoint(worldToSlice, focalPoint, focalPoint);
  if (focalPoint[3] != 0.0)
    {
    focalPoint[0] /= focalPoint[3];
    focalPoint[1] /= focalPoint[3];
    focalPoint[2] /= focalPoint[3];
    }

  // Get the checkerboard spacing and apply the offset fraction
  double checkSpacing[2], checkOffset[2];
  property->GetCheckerboardSpacing(checkSpacing);
  property->GetCheckerboardOffset(checkOffset);
  checkOffset[0] = checkOffset[0]*checkSpacing[0] + focalPoint[0];
  checkOffset[1] = checkOffset[1]*checkSpacing[1] + focalPoint[1];

  // Adjust according to the origin and spacing of the slice data
  double origin[3], spacing[3];
  input->GetSpacing(spacing);
  input->GetOrigin(origin);
  checkOffset[0] = (checkOffset[0] - origin[0])/spacing[0];
  checkOffset[1] = (checkOffset[1] - origin[1])/spacing[1];
  checkSpacing[0] /= spacing[0],
  checkSpacing[1] /= spacing[1];

  // Apply the checkerboard to the data
  int extent[6];
  input->GetExtent(extent);
  unsigned char *data = static_cast<unsigned char *>(
    input->GetScalarPointerForExtent(extent));

  vtkImageMapper3D::CheckerboardRGBA(
    data, extent[1] - extent[0] + 1, extent[3] - extent[2] + 1,
    checkOffset[0], checkOffset[1], checkSpacing[0], checkSpacing[1]);
}

//----------------------------------------------------------------------------
// Compute the vertices of the polygon in the slice coordinate system
void vtkImageResliceMapper::UpdatePolygonCoords(vtkRenderer *ren)
{
  // Get the projection matrix
  double aspect = ren->GetTiledAspectRatio();
  vtkCamera *camera = ren->GetActiveCamera();
  vtkMatrix4x4 *viewMatrix = camera->GetViewTransformMatrix();
  vtkMatrix4x4 *projMatrix = camera->GetProjectionTransformMatrix(
                               aspect, 0, 1);

  // Compute other useful matrices
  double worldToView[16];
  double viewToWorld[16];
  vtkMatrix4x4::Multiply4x4(
    *projMatrix->Element, *viewMatrix->Element, worldToView);
  vtkMatrix4x4::Invert(worldToView, viewToWorld);

  double worldToSlice[16];
  double viewToSlice[16];
  vtkMatrix4x4::Invert(*this->SliceToWorldMatrix->Element, worldToSlice);
  vtkMatrix4x4::Multiply4x4(worldToSlice, viewToWorld, viewToSlice);

  // Get slice plane in world coords by passing null as the matrix
  double plane[4];
  this->GetSlicePlaneInDataCoords(0, plane);

  // Check whether normal is facing towards camera, the "ndop" is
  // the negative of the direction of projection for the camera
  double *ndop = viewMatrix->Element[2];
  if (vtkMath::Dot(ndop, plane) < 0)
    {
    plane[0] = -plane[0];
    plane[1] = -plane[1];
    plane[2] = -plane[2];
    plane[3] = -plane[3];
    }

  // Get the z position of the slice in slice coords
  // (requires plane to be normalized by GetSlicePlaneInDataCoords)
  double z = (plane[2] - 2.0)*plane[3];

  // Generate a tolerance based on the screen pixel size
  double fpoint[4];
  camera->GetFocalPoint(fpoint);
  fpoint[3] = 1.0;
  vtkMatrix4x4::MultiplyPoint(worldToView, fpoint, fpoint);
  fpoint[0] /= fpoint[3];
  fpoint[1] /= fpoint[3];
  fpoint[2] /= fpoint[3];
  fpoint[3] = 1.0;

  double topOfScreen[4], botOfScreen[4];
  fpoint[1] -= 1.0;
  vtkMatrix4x4::MultiplyPoint(viewToWorld, fpoint, topOfScreen);
  fpoint[1] += 2.0;
  vtkMatrix4x4::MultiplyPoint(viewToWorld, fpoint, botOfScreen);

  topOfScreen[0] /= topOfScreen[3];
  topOfScreen[1] /= topOfScreen[3];
  topOfScreen[2] /= topOfScreen[3];
  topOfScreen[3] = 1.0;

  botOfScreen[0] /= botOfScreen[3];
  botOfScreen[1] /= botOfScreen[3];
  botOfScreen[2] /= botOfScreen[3];
  botOfScreen[3] = 1.0;

  // height of view in world coords at focal point
  double viewHeight =
    sqrt(vtkMath::Distance2BetweenPoints(topOfScreen, botOfScreen));

  // height of view in pixels
  int height = ren->GetSize()[1];

  double tol = (height == 0 ? 0.5 : viewHeight*0.5/height); 

  // make the data bounding box (with or without border)
  int border = this->Border;
  double b = (border ? 0.5 : 0.0);
  double bounds[6];
  for (int ii = 0; ii < 3; ii++)
    {
    double c = b*this->DataSpacing[ii];
    int lo = this->DataWholeExtent[2*ii];
    int hi = this->DataWholeExtent[2*ii+1];
    if (border == 0 && lo == hi)
      { // apply tolerance to avoid degeneracy
      c = tol;
      }
    bounds[2*ii]   = lo*this->DataSpacing[ii] + this->DataOrigin[ii] - c;
    bounds[2*ii+1] = hi*this->DataSpacing[ii] + this->DataOrigin[ii] + c;
    }

  // transform the vertices to the slice coord system
  double xpoints[8];
  double ypoints[8];
  double weights[8];
  bool above[8];
  double mat[16];
  vtkMatrix4x4::Multiply4x4(*this->WorldToDataMatrix->Element,
                            *this->SliceToWorldMatrix->Element, mat);
  vtkMatrix4x4::Invert(mat, mat);

  for (int i = 0; i < 8; i++)
    {
    double point[4];
    point[0] = bounds[0 + ((i>>0)&1)];
    point[1] = bounds[2 + ((i>>1)&1)];
    point[2] = bounds[4 + ((i>>2)&1)];
    point[3] = 1.0;
    vtkMatrix4x4::MultiplyPoint(mat, point, point);
    xpoints[i] = point[0]/point[3];
    ypoints[i] = point[1]/point[3];
    weights[i] = point[2]/point[3] - z;
    above[i] = (weights[i] >= 0);
    }

  // go through the edges and find the new points 
  double newxpoints[12];
  double newypoints[12];
  double cx = 0.0;
  double cy = 0.0;
  int n = 0;
  for (int j = 0; j < 12; j++)
    {
    // verts from edges (sorry about this..)
    int i1 = (j & 3) | (((j<<1) ^ (j<<2)) & 4);
    int i2 = (i1 ^ (1 << (j>>2)));

    if (above[i1] ^ above[i2])
      {
      double w1 = weights[i2];
      double w2 = -weights[i1];
      newxpoints[n] = (w1*xpoints[i1] + w2*xpoints[i2])/(w1 + w2);
      newypoints[n] = (w1*ypoints[i1] + w2*ypoints[i2])/(w1 + w2);
      cx += newxpoints[n];
      cy += newypoints[n];
      n++;
      }
    }

  // n should never exceed six
  if (n > 6)
    {
    vtkErrorMacro("UpdateCutPolygon generated more than "
                  "6 points, please report a bug!");
    }

  double *coords = this->Coords;
  this->NCoords = n;

  if (n > 0)
    {
    // centroid 
    cx /= n;
    cy /= n;

    // sort the points to make a convex polygon
    double angles[6];
    for (int k = 0; k < n; k++)
      {
      double x = newxpoints[k];
      double y = newypoints[k];
      double t = atan2(y - cy, x - cx);
      int kk;
      for (kk = 0; kk < k; kk++)
        {
        if (t < angles[kk]) { break; }
        }
      for (int jj = k; jj > kk; --jj)
        {
        int jj3 = jj*3;
        angles[jj] = angles[jj-1];
        coords[jj3] = coords[jj3-3]; 
        coords[jj3+1] = coords[jj3-2]; 
        coords[jj3+2] = coords[jj3-1]; 
        }
      int kk3 = kk*3;
      angles[kk] = t;
      coords[kk3] = x;
      coords[kk3+1] = y;
      coords[kk3+2] = z;
      }
    }
}

//----------------------------------------------------------------------------
// Compute the texture coordinates for the cut polygon
void vtkImageResliceMapper::ComputeTCoords(
  vtkImageData *input, const int extent[6], int ncoords,
  const double *coords, double *tcoords)
{
  // info about the texture, based on the provided extent
  int xdim, ydim;
  int imageSize[2];
  int textureSize[2];

  // compute image size and texture size from extent
  this->ComputeTextureSize(
    extent, xdim, ydim, imageSize, textureSize);

  // now get the info about the image
  double *spacing = input->GetSpacing();
  double *origin = input->GetOrigin();

  // compute the texture coords
  for (int k = 0; k < ncoords; k++)
    {
    int k2 = k*2;
    int k3 = k*3;
    tcoords[k2]   = ((coords[k3]   - origin[0] + 0.5*spacing[0])/
                     (textureSize[0]*spacing[0]));
    tcoords[k2+1] = ((coords[k3+1] - origin[1] + 0.5*spacing[1])/
                     (textureSize[1]*spacing[1]));
    }
}

//----------------------------------------------------------------------------
void vtkImageResliceMapper::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
unsigned long vtkImageResliceMapper::GetMTime()
{
  unsigned long mTime = this->Superclass::GetMTime();

  // Include camera in MTime so that REQUEST_INFORMATION
  // will be called if the camera changes
  if (this->SliceFacesCamera || this->SliceAtFocalPoint)
    {
    vtkRenderer *ren = this->GetCurrentRenderer();
    if (ren)
      {
      vtkCamera *camera = ren->GetActiveCamera();
      unsigned long mTime2 = camera->GetMTime();
      if (mTime2 > mTime)
        {
        mTime = mTime2;
        }
      }
    }

  if (!this->SliceFacesCamera || !this->SliceAtFocalPoint)
    {
    unsigned long sTime = this->SlicePlane->GetMTime();
    if (sTime > mTime)
      {
      mTime = sTime;
      }
    }

  return mTime;
}

//----------------------------------------------------------------------------
double *vtkImageResliceMapper::GetBounds()
{
  // Modify to give just the slice bounds
  if (!this->GetInput())
    {
    vtkMath::UninitializeBounds(this->Bounds);
    return this->Bounds;
    }
  else
    {
    this->UpdateInformation();
    double *spacing = this->DataSpacing;
    double *origin = this->DataOrigin;
    int *extent = this->DataWholeExtent;

    int swapXBounds = (spacing[0] < 0);  // 1 if true, 0 if false
    int swapYBounds = (spacing[1] < 0);  // 1 if true, 0 if false
    int swapZBounds = (spacing[2] < 0);  // 1 if true, 0 if false

    this->Bounds[0] = origin[0] + (extent[0+swapXBounds] * spacing[0]);
    this->Bounds[2] = origin[1] + (extent[2+swapYBounds] * spacing[1]);
    this->Bounds[4] = origin[2] + (extent[4+swapZBounds] * spacing[2]);

    this->Bounds[1] = origin[0] + (extent[1-swapXBounds] * spacing[0]);
    this->Bounds[3] = origin[1] + (extent[3-swapYBounds] * spacing[1]);
    this->Bounds[5] = origin[2] + (extent[5-swapZBounds] * spacing[2]);

    return this->Bounds;
    }
}

//----------------------------------------------------------------------------
void vtkImageResliceMapper::ReportReferences(vtkGarbageCollector* collector)
{
  this->Superclass::ReportReferences(collector);
  // These filters share our input and are therefore involved in a
  // reference loop.
  vtkGarbageCollectorReport(collector, this->ImageReslice, "ImageReslice");
}
