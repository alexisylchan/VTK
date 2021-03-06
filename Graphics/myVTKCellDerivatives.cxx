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
#include <assert.h>


#include <math.h>

#include "vtkMath.h"

//typedef TNT::Array1D<double> ArrayType;
//typedef TNT::Array2D<double> MatrixType;
//typedef JAMA::Eigenvalue<double> SolverType;


vtkStandardNewMacro(myVTKCellDerivatives);
int debug = 0;

myVTKCellDerivatives::myVTKCellDerivatives()
{
  this->VectorMode = VTK_VECTOR_MODE_COMPUTE_GRADIENT;
  this->TensorMode = VTK_TENSOR_MODE_COMPUTE_GRADIENT;
  this->ComputeQ = VTK_NO_Q;
  this->ComputeCurvatureSJ = VTK_NO_CURVATURE_SJ;
  this->ComputeCurvatureHO = VTK_NO_CURVATURE_HO;


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

  vtkDoubleArray *outCurvatureSJ = NULL;
  vtkDoubleArray *outCurvatureSJVector = NULL;
  vtkDoubleArray *outCurvatureHO = NULL;
  vtkDoubleArray *outCurvatureHOVector = NULL;

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
		
		if (this->ComputeCurvatureSJ == VTK_CURVATURE_SJ)  // Technically only 1 scalar allowed, so if ComputeQ is on, this will override ComputeQ
		{
			outCurvatureSJ = vtkDoubleArray::New();
			outCurvatureSJ->SetNumberOfComponents(1);
			outCurvatureSJ->SetNumberOfTuples(numCells);
			outCurvatureSJ->SetName("CurvatureSJ");

			outCurvatureSJVector = vtkDoubleArray::New();
			outCurvatureSJVector->SetNumberOfComponents(3);
			outCurvatureSJVector->SetNumberOfTuples(numCells);
			outCurvatureSJVector->SetName("Acceleration");
		}
		if (this->ComputeCurvatureHO == VTK_CURVATURE_HO)  // Technically only 1 scalar allowed, so if ComputeQ is on, this will override ComputeQ
		{
			outCurvatureHO = vtkDoubleArray::New();
			outCurvatureHO->SetNumberOfComponents(1);
			outCurvatureHO->SetNumberOfTuples(numCells);
			outCurvatureHO->SetName("CurvatureHO");

			outCurvatureHOVector = vtkDoubleArray::New();
			outCurvatureHOVector->SetNumberOfComponents(3);
			outCurvatureHOVector->SetNumberOfTuples(numCells);
			outCurvatureHOVector->SetName("CurvatureHOVector");
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
    double pcoords[3], derivs[9], w[3], *scalars, *vectors, Q[1], curvatureSJ[1], curvatureHO[1];
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
    vtkTensor* velThirdOrderTens = vtkTensor::New();

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

		  //ArrayType sortedRealEigValue = this->getLambda2(realeigvalue);
		  double lambda2[1];
		  lambda2[0] = getMedianEigvalue(realeigvalue, ieigvalue); 
		  outLambda2->SetTuple(cellId,lambda2);

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
		if (this->ComputeCurvatureSJ == VTK_CURVATURE_SJ)
		{
			double acc[3];
			//Jv
			acc[0] = derivs[0]*vectors[0] + derivs[1]*vectors[1] + derivs[2]*vectors[2];
			acc[1] = derivs[3]*vectors[0] + derivs[4]*vectors[1] + derivs[5]*vectors[2];
			acc[2] = derivs[6]*vectors[0] + derivs[7]*vectors[1] + derivs[8]*vectors[2];
			double crossproduct[3];
			
			//v x Jv      or     v x a
			vtkMath::Cross(vectors,acc,crossproduct);
			
			//c =  ( v x Jv )/(|v|^3)
			for (int m = 0; m < 3; m++)
			{
				crossproduct[m] /= pow(vtkMath::Norm(vectors),3);
			}
			
			// c dot v
			curvatureSJ[0] = vtkMath::Dot(crossproduct,vectors);// Get whether velocity is parallel to acceleration
			
			outCurvatureSJ->SetTuple(cellId,curvatureSJ);
			outCurvatureSJVector->SetTuple(cellId,acc);// Debug: Draw acceleration instead of cross product 
			 
		}
		if (this->ComputeCurvatureHO == VTK_CURVATURE_HO)
		{
			//double acc[3];

			////Jv
			//acc[0] = derivs[0]*vectors[0] + derivs[1]*vectors[1] + derivs[2]*vectors[2];
			//acc[1] = derivs[3]*vectors[0] + derivs[4]*vectors[1] + derivs[5]*vectors[2];
			//acc[2] = derivs[6]*vectors[0] + derivs[7]*vectors[1] + derivs[8]*vectors[2];
			//
			////JJv
			//double grad_acc[3]; 
			//grad_acc[0] = derivs[0]*acc[0] + derivs[1]*acc[1] + derivs[2]*acc[2];
			//grad_acc[1] = derivs[3]*acc[0] + derivs[4]*acc[1] + derivs[5]*acc[2];
			//grad_acc[2] = derivs[6]*acc[0] + derivs[7]*acc[1] + derivs[8]*acc[2];
			//
			//double velThirdOrder[27];

			//// T = third-order tensor of second derivatives of v
			////cell->Derivatives(0, pcoords, derivs, 3, velThirdOrder);
			//int thirdOrderCount = 0;
			//for (int k = 0; k<9; k++)
			//{ 
			//	double  thirdOrderTuple[3];
			//	cell->Derivatives(0,pcoords,derivs+k,3,thirdOrderTuple);
			//	for (int m = 0; m < 3; m++)
			//	{
			//		//if (vtkMath::IsNan(thirdOrderTuple[m]) )//|| vtkMath::IsInf(thirdOrderTuple[m]))//Debug: remove extreme values first
			//		//{	
			//		//	velThirdOrder[thirdOrderCount]= 0; thirdOrderCount++;
			//		//}
			//		//else
			//		//{	
			//			velThirdOrder[thirdOrderCount]=thirdOrderTuple[m]; thirdOrderCount++;
			//		/*}*/
			//	}
			//}
			//// Tv
			//double Tv[9];
			//Tv[0] = velThirdOrder[0]*vectors[0] + velThirdOrder[1]*vectors[1] + velThirdOrder[2]*vectors[2];
			//Tv[1] = velThirdOrder[3]*vectors[0] + velThirdOrder[4]*vectors[1] + velThirdOrder[5]*vectors[2];
			//Tv[2] = velThirdOrder[6]*vectors[0] + velThirdOrder[7]*vectors[1] + velThirdOrder[8]*vectors[2];
			//Tv[3] = velThirdOrder[9]*vectors[0] + velThirdOrder[10]*vectors[1] + velThirdOrder[11]*vectors[2];
			//Tv[4] = velThirdOrder[12]*vectors[0] + velThirdOrder[13]*vectors[1] + velThirdOrder[14]*vectors[2];
			//Tv[5] = velThirdOrder[15]*vectors[0] + velThirdOrder[16]*vectors[1] + velThirdOrder[17]*vectors[2];
			//Tv[6] = velThirdOrder[18]*vectors[0] + velThirdOrder[19]*vectors[1] + velThirdOrder[20]*vectors[2];
			//Tv[7] = velThirdOrder[21]*vectors[0] + velThirdOrder[22]*vectors[1] + velThirdOrder[23]*vectors[2];
			//Tv[8] = velThirdOrder[24]*vectors[0] + velThirdOrder[25]*vectors[1] + velThirdOrder[26]*vectors[2];
			//if (debug)
			//	cout<<"Tv "<<" "<<Tv[0]<<" "<<Tv[1]<<" "<<Tv[2]<<" "<<Tv[3]<<" "<<Tv[4]<<" "<<Tv[5]<<" "<<Tv[6]<<" "<<Tv[7]<<" "<<Tv[8]<<endl;
			//
			////Tvv
			//double Tvv[3]; 
			//Tvv[0] = Tv[0]*vectors[0] + Tv[1]*vectors[1] + Tv[2]*vectors[2];
			//Tvv[1] = Tv[3]*vectors[0] + Tv[4]*vectors[1] + Tv[5]*vectors[2];
			//Tvv[2] = Tv[6]*vectors[0] + Tv[7]*vectors[1] + Tv[8]*vectors[2];
			//if (debug)
			//	cout<<"Tvv "<<" "<<Tvv[0]<<" "<<Tvv[1]<<" "<<Tvv[2]<<endl;
			//

			//double thirdDerivOfStreamline[3];
			//// b = JJv + Tvv
			//vtkMath::Add(grad_acc,Tvv,thirdDerivOfStreamline);
			//if (debug)
			//	cout<<"thirdDerivOfStreamline "<<thirdDerivOfStreamline[0]<<" "<<thirdDerivOfStreamline[1]<<" "<<thirdDerivOfStreamline[2]<<endl;
			//
			//
			//double crossproduct[3];	
			//vtkMath::Cross(vectors,thirdDerivOfStreamline,crossproduct);
			//if (debug)
			//	cout<<"crossproduct "<<crossproduct[0]<<" "<<crossproduct[1]<<" "<<crossproduct[2]<<endl;
			//
			//curvatureHO[0] = vtkMath::Norm(crossproduct,3);
			//if (debug)
			//	cout<<"crossproduct mag"<<curvatureHO[0]<<endl;


			////double crossproduct[3];			  
			//////v x Jv      or     v x a
			////vtkMath::Cross(vectors,acc,crossproduct);
			////if (debug)
			////	cout<<"crossproduct "<<crossproduct[0]<<" "<<crossproduct[1]<<" "<<crossproduct[2]<<endl;

			////// torsion = ((v x Jv) dot b) / | (v x Jv)|^2
			////double dotProduct = vtkMath::Dot(crossproduct,thirdDerivOfStreamline);//
			////curvatureHO[0] = dotProduct/pow(vtkMath::Norm(crossproduct),2); // curvatureHO is torsion
			////if (vtkMath::IsNan(curvatureHO[0]) || vtkMath::IsInf(curvatureHO[0]))
			////	curvatureHO[0]=-4.69385*exp(306.0); // set to lowest known value
			////if (debug)
			////{
			////	cout<<"dot"<<dotProduct<<endl;
			////	cout<<"curvatureHO"<<curvatureHO[0]<<endl;
			////}
			//outCurvatureHO->SetTuple(cellId,curvatureHO);
			//outCurvatureHOVector->SetTuple(cellId,thirdDerivOfStreamline); // Draw Jacobian of acceleration
			 
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
  if (outLambda2)
  {
	  outCD->SetScalars(outLambda2);
	  outLambda2->Delete();
  }
  if (outQ && this->ComputeQ == VTK_COMPUTE_Q)
  {
	  outCD->SetScalars(outQ);
	  outQ->Delete();
  }
  if (outCurvatureSJ && this->ComputeCurvatureSJ == VTK_CURVATURE_SJ)
  {
	  outCD->SetScalars(outCurvatureSJ);
	  outCurvatureSJ->Delete();
	  outCD->SetVectors(outCurvatureSJVector);
	  outCurvatureSJVector->Delete();
  }
  if (outCurvatureHO && this->ComputeCurvatureHO == VTK_CURVATURE_HO)
  {
	  outCD->SetScalars(outCurvatureHO);
	  outCurvatureHO->Delete();
	  outCD->SetVectors(outCurvatureHOVector);
	  outCurvatureHOVector->Delete();
  }


  return 1;
}
ArrayType myVTKCellDerivatives::getLambda2(ArrayType realeigvalue)
{ 
	if (realeigvalue.dim()<=1)
		return realeigvalue;
	ArrayType currArray = realeigvalue.subarray(1,realeigvalue.dim()-1);
	ArrayType lessArray(realeigvalue.dim()-1);
	ArrayType greaterArray(realeigvalue.dim()-1);
	int lessArrayInt = 0;
	int greaterArrayInt = 0;
	for (int i =0; i<currArray.dim();i++)
	{
		if (currArray[i] <= realeigvalue[0])
		{
			lessArray[lessArrayInt] = currArray[i];
			lessArrayInt++;
		}
		else
		{
			greaterArray[greaterArrayInt] = currArray[i];
			greaterArrayInt++;
		}
			
	}
	ArrayType returnArray(realeigvalue.dim());
	returnArray.inject(getLambda2(lessArray));
	returnArray.inject(ArrayType(1,realeigvalue[0]));
	returnArray.inject(getLambda2(greaterArray));
	return returnArray;
     //var list less, greater
     //if length(array) ≤ 1
     //    return array  // an array of zero or one elements is already sorted
     //select and remove a pivot value pivot from array
     //for each x in array
     //    if x ≤ pivot then append x to less
     //    else append x to greater
     //return concatenate(quicksort(less), pivot, quicksort(greater))

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

double myVTKCellDerivatives::getMedianEigvalue(ArrayType realeigvalue, ArrayType ieigvalue)
{
		 
		  
		 /* if (realeigvalue.dim() >= 3)
		  {
		   cout<<"Actual "<<realeigvalue[0]<<","<<realeigvalue[1]<<","<<realeigvalue[2]<<endl;
		  }
		  else
		  {
			  cout<<ieigvaluesize<<endl;
		  }*/
		  /*for (int l = 0; l < ieigvalue.dim(); l++)
		  {
			  if (abs(ieigvalue[l]) > 0.1)
				 
		  }*/
		 /* cout<<"Real "<<realeigvalue[0]<<","<<realeigvalue[1]<<","<<realeigvalue[2]<<endl;
		  cout<<"Imag "<<ieigvalue[0]<<","<<ieigvalue[1]<<","<<ieigvalue[2]<<endl;*/
		  for (int i = 0; i < realeigvalue.dim();i++)
		  {
			  for (int j = 0; j < realeigvalue.dim();j++)
			  {
				 for (int k = 0; k < realeigvalue.dim() - 1 - j; k++)
				 {
					 if (realeigvalue[k] > realeigvalue[k+1])
					 {
						 double temp = realeigvalue[k+1];
						 realeigvalue[k+1] = realeigvalue[k];
						 realeigvalue[k] = temp;
					 }
				 }
			  }
		  } 
		  // cout<<"Sorted "<<realeigvalue[0]<<","<<realeigvalue[1]<<","<<realeigvalue[2]<<endl;
		  return realeigvalue[1];
}
