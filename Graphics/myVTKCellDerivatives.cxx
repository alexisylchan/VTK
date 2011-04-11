/*=========================================================================

  Program:   Visualization Toolkit
  Module:    myVTKCellDerivatives.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "myVTKCellDerivatives.h"

#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkDataSet.h"
#include "vtkDoubleArray.h"
#include "vtkGenericCell.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkTensor.h"

#include "jama_eig.h"

#include <math.h>

typedef TNT::Array1D<double> ArrayType;
typedef TNT::Array2D<double> MatrixType;
typedef JAMA::Eigenvalue<double> SolverType;


vtkStandardNewMacro(myVTKCellDerivatives);

myVTKCellDerivatives::myVTKCellDerivatives()
{
  this->VectorMode = VTK_VECTOR_MODE_COMPUTE_GRADIENT;
  this->TensorMode = VTK_TENSOR_MODE_COMPUTE_GRADIENT;
  this->ComputeQ = VTK_NO_Q;

  // by default process active point scalars
  this->SetInputArrayToProcess(0,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::SCALARS);

  // by default process active point vectors
  this->SetInputArrayToProcess(1,0,0,vtkDataObject::FIELD_ASSOCIATION_POINTS,
                               vtkDataSetAttributes::VECTORS);
}

int myVTKCellDerivatives::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkDataSet *input = vtkDataSet::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkDataSet *output = vtkDataSet::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkPointData *pd=input->GetPointData(), *outPD=output->GetPointData();
  vtkCellData *cd=input->GetCellData(), *outCD=output->GetCellData();
  vtkDataArray *inScalars=this->GetInputArrayToProcess(0, inputVector);
  vtkDataArray *inVectors=this->GetInputArrayToProcess(1, inputVector);
  vtkDoubleArray *outGradients=NULL;
  vtkDoubleArray *outVorticity=NULL;
  vtkDoubleArray *outTensors=NULL;
  vtkDoubleArray *outQ=NULL;
  vtkDoubleArray *outLambda2=NULL;
  vtkIdType numCells=input->GetNumberOfCells();
  int computeScalarDerivs=1, computeVectorDerivs=1, computeVorticity=1, subId;

  // Initialize
  vtkDebugMacro(<<"Computing cell derivatives");

  // First, copy the input to the output as a starting point
  output->CopyStructure( input );

  // Check input
  if ( numCells < 1 )
    {
    vtkErrorMacro("No cells to generate derivatives from");
    return 1;
    }

  // Figure out what to compute
  if ( inScalars && this->VectorMode == VTK_VECTOR_MODE_COMPUTE_GRADIENT )
    {
    outGradients = vtkDoubleArray::New();
    outGradients->SetNumberOfComponents(3);
    outGradients->SetNumberOfTuples(numCells);
    outGradients->SetName("ScalarGradient");
		if (this->ComputeQ == VTK_COMPUTE_Q)
		{
			outQ = vtkDoubleArray::New();
			outQ->SetNumberOfComponents(1);
			outQ->SetNumberOfTuples(numCells);
			outQ->SetName("QCriteria");
		}
    }
  else
    {
    computeScalarDerivs = 0;
    }

  if ( inVectors && this->VectorMode == VTK_VECTOR_MODE_COMPUTE_VORTICITY )
    {
    outVorticity = vtkDoubleArray::New();
    outVorticity->SetNumberOfComponents(3);
    outVorticity->SetNumberOfTuples(numCells);
    outVorticity->SetName("Vorticity");
    }
  else
    {
    computeVorticity = 0;
    }

  if (inVectors && ( this->TensorMode == VTK_TENSOR_MODE_COMPUTE_GRADIENT ||
                     this->TensorMode == VTK_TENSOR_MODE_COMPUTE_SPIN|| 
                     this->TensorMode == VTK_TENSOR_MODE_COMPUTE_LAMBDA2||
                     this->TensorMode == VTK_TENSOR_MODE_COMPUTE_STRAIN ))
    {
    outTensors = vtkDoubleArray::New();
    outTensors->SetNumberOfComponents(9);
    outTensors->SetNumberOfTuples(numCells);
    if ( this->TensorMode == VTK_TENSOR_MODE_COMPUTE_STRAIN )
      {
      outTensors->SetName("Strain");
      }
    else if ( this->TensorMode == VTK_TENSOR_MODE_COMPUTE_GRADIENT )
      {
      outTensors->SetName("VectorGradient");
      }
    else if ( this->TensorMode == VTK_TENSOR_MODE_COMPUTE_SPIN )
      {
      outTensors->SetName("Spin");
      }
    else if ( this->TensorMode == VTK_TENSOR_MODE_COMPUTE_LAMBDA2 )
      {
      outTensors->SetName("Lambda2Tensors");
	  outLambda2 = vtkDoubleArray::New();
	  outLambda2->SetNumberOfComponents(1);
	  outLambda2->SetNumberOfTuples(numCells);
	  outLambda2->SetName("Lambda2");
      }
    }
  else
    {
    computeVectorDerivs = 0;
    }

  // If just passing data forget the loop
  if ( computeScalarDerivs || computeVectorDerivs || computeVorticity )
    {
    double pcoords[3], derivs[9], w[3], *scalars, *vectors, Q[1];
    vtkGenericCell *cell = vtkGenericCell::New();
    vtkIdType cellId;
    vtkDoubleArray *cellScalars=vtkDoubleArray::New();
    if ( computeScalarDerivs )
      {
      cellScalars->SetNumberOfComponents(inScalars->GetNumberOfComponents());
      cellScalars->Allocate(cellScalars->GetNumberOfComponents()*VTK_CELL_SIZE);
      cellScalars->SetName("Scalars");
      }
    vtkDoubleArray *cellVectors=vtkDoubleArray::New(); 
    cellVectors->SetNumberOfComponents(3);
    cellVectors->Allocate(3*VTK_CELL_SIZE);
    cellVectors->SetName("Vectors");
    vtkTensor* tens = vtkTensor::New();

    // Loop over all cells computing derivatives
    vtkIdType progressInterval = numCells/20 + 1;
    for (cellId=0; cellId < numCells; cellId++)
      {
      if ( ! (cellId % progressInterval) ) 
        {
        vtkDebugMacro(<<"Computing cell #" << cellId);
        this->UpdateProgress (static_cast<double>(cellId)/numCells);
        }

      input->GetCell(cellId, cell);
      subId = cell->GetParametricCenter(pcoords);
      
      if ( computeScalarDerivs )
        {
        inScalars->GetTuples(cell->PointIds, cellScalars);
        scalars = cellScalars->GetPointer(0);
		//((d(vx)/dx),(d(vx)/dy),(d(vx)/dz), (d(vy)/dx),(d(vy)/dy), (d(vy)/dz), (d(vz)/dx),(d(vz)/dy),(d(vz)/dz)). 
        cell->Derivatives(subId, pcoords, scalars, 1, derivs); 
        outGradients->SetTuple(cellId, derivs);
        }

      if ( computeVectorDerivs || computeVorticity )
        {
        inVectors->GetTuples(cell->PointIds, cellVectors);
        vectors = cellVectors->GetPointer(0);
        cell->Derivatives(0, pcoords, vectors, 3, derivs);

        // Insert appropriate tensor
        if ( this->TensorMode == VTK_TENSOR_MODE_COMPUTE_GRADIENT)
          {
          tens->SetComponent(0,0, derivs[0]);
          tens->SetComponent(0,1, derivs[1]);
          tens->SetComponent(0,2, derivs[2]);
          tens->SetComponent(1,0, derivs[3]);
          tens->SetComponent(1,1, derivs[4]);
          tens->SetComponent(1,2, derivs[5]);
          tens->SetComponent(2,0, derivs[6]);
          tens->SetComponent(2,1, derivs[7]);
          tens->SetComponent(2,2, derivs[8]);
          
          outTensors->InsertTuple(cellId, tens->T);
          }
        else if(this->TensorMode == VTK_TENSOR_MODE_COMPUTE_STRAIN)
          {
          tens->SetComponent(0,0, derivs[0]);
          tens->SetComponent(0,1, 0.5*(derivs[1]+derivs[3]));
          tens->SetComponent(0,2, 0.5*(derivs[2]+derivs[6]));
          tens->SetComponent(1,0, 0.5*(derivs[1]+derivs[3]));
          tens->SetComponent(1,1, derivs[4]);
          tens->SetComponent(1,2, 0.5*(derivs[5]+derivs[7]));
          tens->SetComponent(2,0, 0.5*(derivs[2]+derivs[6]));
          tens->SetComponent(2,1, 0.5*(derivs[5]+derivs[7]));
          tens->SetComponent(2,2, derivs[8]);
          
          outTensors->InsertTuple(cellId, tens->T);
          }
        else if(this->TensorMode == VTK_TENSOR_MODE_COMPUTE_SPIN)
          {
          tens->SetComponent(0,0, derivs[0]);
          tens->SetComponent(0,1, 0.5*(derivs[1]-derivs[3]));
          tens->SetComponent(0,2, 0.5*(derivs[2]-derivs[6]));
          tens->SetComponent(1,0, 0.5*(derivs[1]-derivs[3]));
          tens->SetComponent(1,1, derivs[4]);
          tens->SetComponent(1,2, 0.5*(derivs[5]-derivs[7]));
          tens->SetComponent(2,0, 0.5*(derivs[2]-derivs[6]));
          tens->SetComponent(2,1, 0.5*(derivs[5]-derivs[7]));
          tens->SetComponent(2,2, derivs[8]);
          
          outTensors->InsertTuple(cellId, tens->T);
          }
		//else if (this->TensorMode == VTK_TENSOR_MODE_COMPUTE_LAMBDA2)
		//{
		//	Q=0.5*((-derivs[0]*derivs[0]-derivs[4]*derivs[4]-derivs[8]*derivs[8])-2*(derivs[1]*derivs[3]+derivs[2]*derivs[6]+derivs[5]*derivs[7]));
		//}
		//https://cei-gate.ceintl.com/bugs-ext/show_bug.cgi?id=3875
		//((d(vx)/dx),(d(vx)/dy),(d(vx)/dz), (d(vy)/dx),(d(vy)/dy), (d(vy)/dz), (d(vz)/dx),(d(vz)/dy),(d(vz)/dz)). 
		// Q = .5*(-GVxx^2 - GVyy^2 - GVzz^2 - 2(GVxy GVyx + GVxz GVzx + GVyz GVzy))

        else if(this->TensorMode == VTK_TENSOR_MODE_COMPUTE_LAMBDA2)// This is not actually Lambda2, but the tensor used for calculating lambda2
          {
		  MatrixType eigmatrix(3,3);

          tens->SetComponent(0,0, derivs[0]*derivs[0] +  0.25*derivs[1]*derivs[1] + 0.25*derivs[3]*derivs[3] + 0.25*derivs[2]*derivs[2] + 0.25*derivs[6]*derivs[6]);
		  eigmatrix[0][0] = tens->GetComponent(0,0);
		  tens->SetComponent(0,1, 0.5*derivs[0]*derivs[1] + 0.5*derivs[4]*derivs[1] + 0.25*derivs[2]*derivs[5] + 0.25*derivs[6]*derivs[7]);
          eigmatrix[0][1] = tens->GetComponent(0,1);
		  tens->SetComponent(0,2, 0.5*derivs[0]*derivs[2] + 0.25*derivs[1]*derivs[5] + 0.25*derivs[3]*derivs[7] + 0.5*derivs[8]*derivs[1]);
          eigmatrix[0][2] = tens->GetComponent(0,2);
		  tens->SetComponent(1,0, 0.5*derivs[0]*derivs[1] + 0.5*derivs[4]*derivs[1] + 0.25*derivs[5]*derivs[2] + 0.25*derivs[7]*derivs[6]);
          eigmatrix[1][0] = tens->GetComponent(1,0);
		  tens->SetComponent(1,1, 0.25*derivs[1]*derivs[1] + 0.25*derivs[3]*derivs[3] + 1*derivs[4]*derivs[4] + 0.25*derivs[5]*derivs[5] + 0.25*derivs[7]*derivs[7]);
          eigmatrix[1][1] = tens->GetComponent(1,1);
		  tens->SetComponent(1,2, 0.25*derivs[1]*derivs[2] + 0.25*derivs[3]*derivs[6] + 0.5*derivs[4]*derivs[5] + 0.5*derivs[8]*derivs[5]);
          eigmatrix[1][2] = tens->GetComponent(1,2);
		  tens->SetComponent(2,0, 0.5*derivs[1]*derivs[2] + 0.25*derivs[5]*derivs[1] + 0.25*derivs[7]*derivs[5] + 0.5*derivs[8]*derivs[2]);
          eigmatrix[2][0] = tens->GetComponent(2,0);
		  tens->SetComponent(2,1, 0.25*derivs[2]*derivs[1] + 0.25*derivs[6]*derivs[5] + 0.5*derivs[4]*derivs[5] + 0.5*derivs[8]*derivs[5]);
          eigmatrix[2][1] = tens->GetComponent(2,1);
		  tens->SetComponent(2,2, 0.25*derivs[2]*derivs[2] + 0.25*derivs[6]*derivs[6] + 0.25*derivs[5]*derivs[5] + 0.25*derivs[7]*derivs[7] + 1*derivs[8]*derivs[8]);
          eigmatrix[2][2] = tens->GetComponent(2,2);
		  
		  //Compute Eigenvalue?
		  ArrayType realeigvalue,ieigvalue;
		  SolverType eigsolver(eigmatrix);
		  eigsolver.getRealEigenvalues(realeigvalue);
		  eigsolver.getImagEigenvalues(ieigvalue);
		  double pot_lambda2[1];
		  for (int i = 0; i< realeigvalue.dim();i++)
		  {
			  if (i == 0)
				  pot_lambda2[0] = realeigvalue[i];
			  else if (realeigvalue[i] > pot_lambda2[0])
				  pot_lambda2[0] = realeigvalue[i];
		  }
		  outLambda2->SetTuple(cellId,pot_lambda2);

		  // End Compute Eigenvalue
          outTensors->InsertTuple(cellId, tens->T);
          }

        if ( computeVorticity )
          {
          w[0] = derivs[7] - derivs[5];
          w[1] = derivs[2] - derivs[6];
          w[2] = derivs[3] - derivs[1];
          outVorticity->SetTuple(cellId, w);
          }
		if (this->ComputeQ == VTK_COMPUTE_Q )
		{
			Q[0]=0.5*((-derivs[0]*derivs[0]-derivs[4]*derivs[4]-derivs[8]*derivs[8])-2*(derivs[1]*derivs[3]+derivs[2]*derivs[6]+derivs[5]*derivs[7]));
			outQ->SetTuple(cellId,Q);
		}
        }
      }//for all cells

    cell->Delete();
    cellScalars->Delete();
    cellVectors->Delete();
    tens->Delete();
    }//if something to compute

  // Pass appropriate data through to output
  outPD->PassData(pd);
  outCD->PassData(cd);
  if (outGradients)
    {
    outCD->SetVectors(outGradients);
    outGradients->Delete();
    }
  if (outVorticity)
    {
    outCD->SetVectors(outVorticity);
    outVorticity->Delete();
    }
  if (outTensors)
    {
    outCD->SetTensors(outTensors);
    outTensors->Delete();
    }
  if (outQ && this->ComputeQ == VTK_COMPUTE_Q)
  {
	  outCD->SetScalars(outQ);
	  outQ->Delete();
  }
  if (outLambda2)
  {
	  outCD->SetScalars(outLambda2);
	  outLambda2->Delete();
  }

  return 1;
}

const char *myVTKCellDerivatives::GetVectorModeAsString(void)
{
  if ( this->VectorMode == VTK_VECTOR_MODE_PASS_VECTORS )
    {
    return "PassVectors";
    }
  else if ( this->VectorMode == VTK_VECTOR_MODE_COMPUTE_GRADIENT )
    {
    return "ComputeGradient";
    }
  else //VTK_VECTOR_MODE_COMPUTE_VORTICITY
    {
    return "ComputeVorticity";
    }
}

const char *myVTKCellDerivatives::GetTensorModeAsString(void)
{
  if ( this->TensorMode == VTK_TENSOR_MODE_PASS_TENSORS )
    {
    return "PassTensors";
    }
  else if ( this->TensorMode == VTK_TENSOR_MODE_COMPUTE_GRADIENT )
    {
    return "ComputeGradient";
    }
  else if (this->TensorMode == VTK_TENSOR_MODE_COMPUTE_STRAIN)
    {
    return "ComputeStrain";
    }
  else if (this->TensorMode == VTK_TENSOR_MODE_COMPUTE_SPIN)
    {
    return "ComputeSpin";
    }
  else if (this->TensorMode == VTK_TENSOR_MODE_COMPUTE_LAMBDA2)
  {
	  return "ComputeLambda2";
  }
}

void myVTKCellDerivatives::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Vector Mode: " << this->GetVectorModeAsString() 
     << endl;

  os << indent << "Tensor Mode: " << this->GetTensorModeAsString() 
     << endl;
}

