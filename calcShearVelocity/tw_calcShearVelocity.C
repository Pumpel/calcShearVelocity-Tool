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

    // Adding the options for commmand line tools in OpenFOAM
    timeSelector::addOptions(true, true);
    Info<< "Times found:" << runTime.times() << endl;

    // Creating a list with all time steps selected via the command line option
    instantList timeDirs = timeSelector::select0(runTime, args);
    Info<< "Times selected:" << timeDirs << endl;
    Info<< "\nEnd\n" << endl;

    forAll(timeDirs, timei)
    {
    	Info << "Time directory: " << timeDirs[timei].value() << endl;
    	// Setting the current time step that will be worked on
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
		volTensorField shearVelocity
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
			dimensionedTensor
			(
				"shearVelocity",
				dimensionSet(0, 0, -1, 0, 0, 0, 0),
				tensor (0,0,0,0,0,0,0,0,0)
			)
		);

		// Jacobian and the transposed Jacobian matrix of the velocity field
		volTensorField J(fvc::grad(U));
		volTensorField J_T(fvc::grad(U));

		// Symmetric part of the Jacobian matrix J, or rate-of-strain tensor
		volTensorField E (0.5*(J + J_T));

		//Decomposition of the symmetric part E
		  // E = Tr(E) + E - Tr(E) d_ij
		  // Where E - Tr(E) is rate-of-shear, the deviatoric part (of the rate-of-strain tensor,
		  // which is in turn the Jacobian matrix of the velocity field
		  // and Tr(E) d_ij the rate-of-expansion tensor, the hydrostatic part (of the rate-of-strain tensor)
		// Here, the new term shear velocity tensor is used. It was formely named shear rate tensor or rate-of-shear
		// tensor
		shearVelocity = (E - 1/3*tr(2*E)*tensor::one);
		shearVelocity.write();


    }

    return 0;
}


// ************************************************************************* //
