/*=========================================================================

  Program:   Visualization Toolkit
  Module:    myVTKCellDerivatives.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME myVTKCellDerivatives - compute derivatives of scalars and vectors
// .SECTION Description
// myVTKCellDerivatives is a filter that computes derivatives of scalars
// and vectors at the center of cells. You can choose to generate
// different output including the scalar gradient (a vector), computed
// tensor vorticity (a vector), gradient of input vectors (a tensor),
// and strain matrix of the input vectors (a tensor); or you may
// choose to pass data through to the output.
//
// Note that it is assumed that on input scalars and vector point data
// is available, which are then used to generate cell vectors and tensors.
// (The interpolation functions of the cells are used to compute the
// derivatives which is why point data is required.)

// .SECTION Caveats
// The computed derivatives are cell attribute data; you can convert them to
// point attribute data by using the vtkCellDataToPointData filter.
// Note that, due to the interpolation function used (obtained using 
// 1/r**2 normalized sum), the derivatives calculated for polygons
// with more than 4 vertices are inaccurate in most cases.
//
// The point data is passed through the filter to the output.

// .SECTION See Also
// vtkVectorNorm

#ifndef __myVTKCellDerivatives_h
#define __myVTKCellDerivatives_h

#include "vtkDataSetAlgorithm.h"

#define VTK_VECTOR_MODE_PASS_VECTORS      0
#define VTK_VECTOR_MODE_COMPUTE_GRADIENT  1
#define VTK_VECTOR_MODE_COMPUTE_VORTICITY 2

#define VTK_TENSOR_MODE_PASS_TENSORS     0
#define VTK_TENSOR_MODE_COMPUTE_GRADIENT 1
#define VTK_TENSOR_MODE_COMPUTE_STRAIN   2
#define VTK_TENSOR_MODE_COMPUTE_SPIN    3
#define VTK_TENSOR_MODE_COMPUTE_LAMBDA2       4

#define VTK_NO_Q 0
#define VTK_COMPUTE_Q 1
 
#include "jama_eig.h"
typedef TNT::Array1D<double> ArrayType;
typedef TNT::Array2D<double> MatrixType;
typedef JAMA::Eigenvalue<double> SolverType;

class  VTK_GRAPHICS_EXPORT myVTKCellDerivatives : public vtkDataSetAlgorithm 
{
public:
  vtkTypeMacro(myVTKCellDerivatives,vtkDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Construct to compute the gradient of the scalars and vectors.
  static myVTKCellDerivatives *New();

  // Description:
  // Control how the filter works to generate vector cell data. You
  // can choose to pass the input cell vectors, compute the gradient
  // of the input scalars, or extract the vorticity of the computed
  // vector gradient tensor. By default (VectorModeToComputeGradient),
  // the filter will take the gradient of the input scalar data.
  vtkSetMacro(VectorMode,int);
  vtkGetMacro(VectorMode,int);
  //Description:
  //Set whether to compute the Q Criteria when VectorMode is set to Compute Gradient.
  vtkSetMacro(ComputeQ,int);
  vtkGetMacro(ComputeQ,int);
  void SetVectorModeToPassVectors() 
    {this->SetVectorMode(VTK_VECTOR_MODE_PASS_VECTORS);};
  void SetVectorModeToComputeGradient() 
    {this->SetVectorMode(VTK_VECTOR_MODE_COMPUTE_GRADIENT);};
  void SetVectorModeToComputeVorticity() 
    {this->SetVectorMode(VTK_VECTOR_MODE_COMPUTE_VORTICITY);};
  void SetComputeQ()
  { this->SetComputeQ(VTK_COMPUTE_Q);};
  void UnsetComputeQ()
  { this->SetComputeQ(VTK_NO_Q);};

  const char *GetVectorModeAsString();

  // Description:
  // Control how the filter works to generate tensor cell data. You can
  // choose to pass the input cell tensors, compute the gradient of
  // the input vectors, or compute the strain tensor of the vector gradient
  // tensor. By default (TensorModeToComputeGradient), the filter will
  // take the gradient of the vector data to construct a tensor.
  vtkSetMacro(TensorMode,int);
  vtkGetMacro(TensorMode,int);
  void SetTensorModeToPassTensors() 
    {this->SetTensorMode(VTK_TENSOR_MODE_PASS_TENSORS);};
  void SetTensorModeToComputeGradient() 
    {this->SetTensorMode(VTK_TENSOR_MODE_COMPUTE_GRADIENT);};
  void SetTensorModeToComputeStrain() 
    {this->SetTensorMode(VTK_TENSOR_MODE_COMPUTE_STRAIN);};
  void SetTensorModeToComputeSpin() 
    {this->SetTensorMode(VTK_TENSOR_MODE_COMPUTE_SPIN);};
  void SetTensorModeToComputeLambda2() 
    {this->SetTensorMode(VTK_TENSOR_MODE_COMPUTE_LAMBDA2);};
  const char *GetTensorModeAsString();

protected:
  myVTKCellDerivatives();
  ~myVTKCellDerivatives() {};
  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  int VectorMode;
  int TensorMode;
  int ComputeQ;
private:
  myVTKCellDerivatives(const myVTKCellDerivatives&);  // Not implemented.
  void operator=(const myVTKCellDerivatives&);  // Not implemented
  ArrayType getLambda2(ArrayType realeigvalue); 
};

#endif
