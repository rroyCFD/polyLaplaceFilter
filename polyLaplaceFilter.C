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

Foam::polyLaplaceFilter::polyLaplaceFilter(const fvMesh& mesh, scalar widthCoeff)
:
    LESfilter(mesh),
    widthCoeff_(widthCoeff),
    coeff_
    (
        IOobject
        (
            "polyLaplaceFilterCoeff",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("zero", dimLength*dimLength, 0),
        calculatedFvPatchScalarField::typeName
    ),

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
    coeff_.ref() = pow(mesh.V(), 2.0/3.0)/widthCoeff_;


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

    deltaSquared_.write();
}


Foam::polyLaplaceFilter::polyLaplaceFilter(const fvMesh& mesh, const dictionary& bd)
:
    LESfilter(mesh),
    widthCoeff_
    (
        readScalar(bd.optionalSubDict(type() + "Coeffs").lookup("widthCoeff"))
    ),
    coeff_
    (
        IOobject
        (
            "polyLaplaceFilterCoeff",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("zero", dimLength*dimLength, 0),
        calculatedFvPatchScalarField::typeName
    ),

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
    coeff_.ref() = pow(mesh.V(), 2.0/3.0)/widthCoeff_;


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

    deltaSquared_.write();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyLaplaceFilter::read(const dictionary& bd)
{
    bd.optionalSubDict(type() + "Coeffs").lookup("widthCoeff") >> widthCoeff_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::polyLaplaceFilter::operator()
(
    const tmp<volScalarField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volScalarField> filteredField =
        unFilteredField() + fvc::laplacian(coeff_, unFilteredField());

    unFilteredField.clear();

    return filteredField;
}


Foam::tmp<Foam::volVectorField> Foam::polyLaplaceFilter::operator()
(
    const tmp<volVectorField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volVectorField> filteredField =
        unFilteredField() + fvc::laplacian(coeff_, unFilteredField());

    unFilteredField.clear();

    return filteredField;
}


Foam::tmp<Foam::volSymmTensorField> Foam::polyLaplaceFilter::operator()
(
    const tmp<volSymmTensorField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volSymmTensorField> filteredField =
        unFilteredField() + fvc::laplacian(coeff_, unFilteredField());

    unFilteredField.clear();

    return filteredField;
}


Foam::tmp<Foam::volTensorField> Foam::polyLaplaceFilter::operator()
(
    const tmp<volTensorField>& unFilteredField
) const
{
    correctBoundaryConditions(unFilteredField);

    tmp<volTensorField> filteredField =
        unFilteredField() + fvc::laplacian(coeff_, unFilteredField());

    unFilteredField.clear();

    return filteredField;
}


// ************************************************************************* //
