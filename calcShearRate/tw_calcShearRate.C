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
    tw_calcShearRate

    The shear rate is calculated as
    	\dot{gamma} = 0.5 (grad(U) + (grad(U))^T)
    where ()^T is the transposed.

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

		Info << "Creating the field shearRate" << endl;
		volSymmTensorField shearRate
		(
			IOobject
			(
				"shearRate",
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionedSymmTensor
			(
				"shearRate",
				dimensionSet(0, 0, -1, 0, 0, 0, 0),
				symmTensor (0,0,0,0,0,0)
			)
		);

		volTensorField myU(fvc::grad(U));
		volTensorField myU_T(myU.T());
		shearRate = symm(myU + myU_T);

		// The diagonal entries of the 2nd rank tensor 'shearRate'
		// must be 0.0
		forAll(shearRate, celli) {
			shearRate[celli].xx() = 0.0;
			shearRate[celli].yy() = 0.0;
			shearRate[celli].zz() = 0.0;
		}

		shearRate.write();
    }

    return 0;
}


// ************************************************************************* //
