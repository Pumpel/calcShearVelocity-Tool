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



		volTensorField myU(fvc::grad(U));
		volTensorField myU_T(myU.T());

//		volSymmTensorField shearVelocity_D
//		(
//			IOobject
//			(
//				"shearVelocity_D",
//				runTime.timeName(),
//				mesh,
//				IOobject::NO_READ,
//				IOobject::AUTO_WRITE
//			),
//			mesh,
//			dimensionedSymmTensor
//			(
//				"shearVelocity_D",
//				dimensionSet(0, 0, -1, 0, 0, 0, 0),
//				symmTensor (0,0,0,0,0,0)
//			)
//		);
//		// Verzerrungstensor (symmetrischen und antisymmetrischen kombiniert)
//		shearVelocity_D = symm(0.5 * (myU + myU_T) );
//		// Obtaining the symmetric part of a tensor
//		// http://mathworld.wolfram.com/SymmetricPart.html
//		shearVelocity_D = (0.5 * (shearVelocity + shearVelocity.T()));
//		// The diagonal entries of the 2nd rank tensor 'shearVelocity'
//		// must be 0.0
//		forAll(shearVelocity, celli) {
//			shearVelocity_D[celli].xx() = 0.0;
//			shearVelocity_D[celli].yy() = 0.0;
//			shearVelocity_D[celli].zz() = 0.0;
//		}
//		shearVelocity_D.write();


		// Shortcut to calculate the shear velocity from the off diagonal
		// elements
		shearVelocity = symm(myU + myU_T);

		// The diagonal entries of the 2nd rank tensor 'shearVelocity'
		// must be 0.0
		forAll(shearVelocity, celli) {
			shearVelocity[celli].xx() = 0.0;
			shearVelocity[celli].yy() = 0.0;
			shearVelocity[celli].zz() = 0.0;
		}

		shearVelocity.write();
    }

    return 0;
}


// ************************************************************************* //
