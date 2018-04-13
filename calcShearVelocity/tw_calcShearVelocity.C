/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    tw_calcShearVelocity

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "transformGeometricField.H"

#include "argList.H"
//#include "Time.H"
#include "timeSelector.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{


    argList::noParallel();
    timeSelector::addOptions();

    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"


    timeSelector::addOptions(true, true);

    Info<< "Times found:" << runTime.times() << endl;

    instantList timeDirs = timeSelector::select0(runTime, args);

    Info<< "Times selected:" << timeDirs << endl;
    Info<< "\nEnd\n" << endl;


    forAll(timeDirs, timei)
    {
    	Info << "Time directory: " << timeDirs[timei].value() << endl;
    	runTime.setTime(timeDirs[timei], timei);
    	mesh.readUpdate();

		Info<< "Reading field U\n" << endl;
		volVectorField U
		(
			IOobject
			(
				"U",
				runTime.timeName(),
				mesh,
				IOobject::MUST_READ,
				IOobject::AUTO_WRITE
			),
			mesh
		);

		Info << "Creating the field shearVelocity" << endl;
		volSymmTensorField shearVelocity
		(
			IOobject
			(
				"shearVelocity",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionedSymmTensor
			(
				"shearVelocity",
				dimensionSet(0, 0, -1, 0, 0, 0, 0),
				symmTensor (0,0,0,0,0,0)
			)
		);

		Info << "Creating the field shearVelocity_E" << endl;
		volScalarField shearVelocity_E
		(
			IOobject
			(
				"shearVelocity_E",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionedScalar
			(
				"shearVelocity_E",
				dimensionSet(0, 0, -1, 0, 0, 0, 0),
				scalar (0)
			)
		);

		Info << "Creating the field E" << endl;
		volTensorField E
		(
			IOobject
			(
				"E",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionedTensor
			(
				"E",
				dimensionSet(0, 0, -1, 0, 0, 0, 0),
				tensor (0,0,0,0,0,0,0,0,0)
			)
		);

		Info << "Creating the field shearVelocity_E_mTr" << endl;
		volTensorField shearVelocity_E_mTr
		(
			IOobject
			(
				"shearVelocity_E_mTr",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionedTensor
			(
				"shearVelocity_E_mTr",
				dimensionSet(0, 0, -1, 0, 0, 0, 0),
				tensor (0,0,0,0,0,0,0,0,0)
			)
		);


		Info << "Creating the field shearVelocity_E_dev" << endl;
		volTensorField shearVelocity_E_dev
		(
			IOobject
			(
				"shearVelocity_E_dev",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionedTensor
			(
				"shearVelocity_E_dev",
				dimensionSet(0, 0, -1, 0, 0, 0, 0),
				tensor (0,0,0,0,0,0,0,0,0)
			)
		);

		// Jacobian and the transposed Jacobian matrix of the velocity field
		volTensorField J(fvc::grad(U));
		volTensorField J_T(fvc::grad(U));

		//********************************************************************************
		// Three different calculation strategies are tested for the rate-of-shear
		// It seems, that the OpenFOAM function symm() doesn't lead to the same (correct)
		// result as the other calculation methods
		//********************************************************************************


		// 1)
		// Testing the symm() function on the Jacobian of the velocity field

		// This calculates the symmetric part, i.e. the rate-of-strain,
		// of the Jacobian of the velocity field
		shearVelocity = symm(J);

		// dev() calculates the deviatoric part, i.e. the rate-of-shear part
		// of the rate-of-strain tensor
		shearVelocity = dev(shearVelocity);

		shearVelocity.write();

		// Doesnt't yield correct results.



		// 2)
		// Testing the result, if all equations are 'handwritten'

		// Symmetric part of the Jacobian matrix J, or rate-of-strain tensor
		E = 0.5*(J + J_T);

		//Decomposition of the symmetric part E
		  // E = Tr(E) + E - Tr(E) d_ij
		  // Where E - Tr(E) is rate-of-shear, the deviatoric part (of the rate-of-strain tensor,
		  // which is in turn the Jacobian matrix of the velocity field
		  // and Tr(E) d_ij the rate-of-expansion tensor, the hydrostatic part (of the rate-of-strain tensor)
		// The two expressions in 2) and 3) should be equal.
		shearVelocity_E_mTr = (E - 1/3*tr(2*E)*tensor::one);

		shearVelocity_E_mTr.write();

		// Yields correct results.


		// 3)
		// Testing the function for the deviatoric part of a tensor

		// dev() calculates the deviatoric part, i.e. the rate-of-shear part
		// of the rate-of-strain tensor
		shearVelocity_E_dev = dev(E);

		shearVelocity_E_dev.write();

		// Yields correct results.


		// 4)
		// Testing the calcluation method for the shear rate according to the lecture notes
		// Hinch, J., Lecture 2: Constitutive Relations

		// The double inner scalar product of the strain rate tensor
		shearVelocity_E = sqrt(2 * E && E);

		shearVelocity_E.write();

		// Doesn't yield correct results.
    }

    return 0;
}


// ************************************************************************* //
