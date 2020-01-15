/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "polyLaplaceFilter.H"
#include "addToRunTimeSelectionTable.H"
#include "calculatedFvPatchFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyLaplaceFilter, 0);
    addToRunTimeSelectionTable(LESfilter, polyLaplaceFilter, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyLaplaceFilter::polyLaplaceFilter(const fvMesh& mesh, scalar d1, scalar d2)
:
    LESfilter(mesh),
    d1_(d1),
    d2_(d2),
    deltaSquared_
    (
        IOobject
        (
            "deltaSquared",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimLength*dimLength, 0)
    )
{
    deltaSquared_.primitiveFieldRef() = (mesh.delta() & mesh.delta());

    forAll(mesh.boundaryMesh(), patchI)
    {
        fvsPatchScalarField& deltaSquaredB =
            deltaSquared_.boundaryFieldRef()[patchI];

        const fvPatch& cPatch = deltaSquaredB.patch();

        if(cPatch.type() != "empty")
        {
            deltaSquaredB = cPatch.delta() & cPatch.delta();
        }
    }

    // deltaSquared_.write();
}


Foam::polyLaplaceFilter::polyLaplaceFilter(const fvMesh& mesh, const dictionary& bd)
:
    LESfilter(mesh),
    d1_(readScalar(bd.optionalSubDict(type() + "Coeffs").lookup("d1"))),
    d2_(readScalar(bd.optionalSubDict(type() + "Coeffs").lookup("d2"))),
    deltaSquared_
    (
        IOobject
        (
            "deltaSquared",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimLength*dimLength, 0)
    )
{
    deltaSquared_.primitiveFieldRef() = (mesh.delta() & mesh.delta());

    forAll(mesh.boundaryMesh(), patchI)
    {
        fvsPatchScalarField& deltaSquaredB =
            deltaSquared_.boundaryFieldRef()[patchI];

        const fvPatch& cPatch = deltaSquaredB.patch();

        if(cPatch.type() != "empty")
        {
            deltaSquaredB = cPatch.delta() & cPatch.delta();
        }
    }

    // deltaSquared_.write();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyLaplaceFilter::read(const dictionary& bd)
{
    bd.optionalSubDict(type() + "Coeffs").lookup("d1") >> d1_;
    bd.optionalSubDict(type() + "Coeffs").lookup("d2") >> d2_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::polyLaplaceFilter::operator()
(
    const tmp<volScalarField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volScalarField> tLapField = fvc::laplacian(deltaSquared_, unFilteredField());

    tmp<volScalarField> filteredField =
        unFilteredField() + d1_*tLapField + d2_*fvc::laplacian(deltaSquared_, tLapField());

    tLapField.clear();
    unFilteredField.clear();

    return filteredField;
}


Foam::tmp<Foam::volVectorField> Foam::polyLaplaceFilter::operator()
(
    const tmp<volVectorField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volVectorField> tLapField = fvc::laplacian(deltaSquared_, unFilteredField());

    tmp<volVectorField> filteredField =
        unFilteredField() + d1_*tLapField + d2_*fvc::laplacian(deltaSquared_, tLapField());

    tLapField.clear();
    unFilteredField.clear();

    return filteredField;
}


Foam::tmp<Foam::volSymmTensorField> Foam::polyLaplaceFilter::operator()
(
    const tmp<volSymmTensorField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volSymmTensorField> tLapField = fvc::laplacian(deltaSquared_, unFilteredField());

    tmp<volSymmTensorField> filteredField =
        unFilteredField() + d1_*tLapField + d2_*fvc::laplacian(deltaSquared_, tLapField());

    tLapField.clear();
    unFilteredField.clear();

    return filteredField;
}


Foam::tmp<Foam::volTensorField> Foam::polyLaplaceFilter::operator()
(
    const tmp<volTensorField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volTensorField> tLapField = fvc::laplacian(deltaSquared_, unFilteredField());

    tmp<volTensorField> filteredField =
        unFilteredField() + d1_*tLapField + d2_*fvc::laplacian(deltaSquared_, tLapField());

    tLapField.clear();
    unFilteredField.clear();

    return filteredField;
}

// ************************************************************************* //
