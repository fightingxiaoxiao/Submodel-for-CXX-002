/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "limitScalarByAlpha.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(limitScalarByAlpha, 0);
    addToRunTimeSelectionTable(option, limitScalarByAlpha, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::limitScalarByAlpha::limitScalarByAlpha
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    scalarName_(coeffs_.get<word>("scalarName")),
    max_(coeffs_.get<scalar>("max")),
    phaseName_(coeffs_.get<word>("phaseName"))
{
    fieldNames_.setSize(1, scalarName_);
    applied_.setSize(1, false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::limitScalarByAlpha::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.readEntry("max", max_);
        coeffs_.readEntry("scalarName", scalarName_);
        coeffs_.readEntry("phaseName", phaseName_);
        return true;
    }

    return false;
}


void Foam::fv::limitScalarByAlpha::correct(volScalarField& S)
{
    const scalar maxScalar = max_;

    volScalarField alpha_ = mesh_.lookupObject<volScalarField>("alpha.air");

    label nLimit = 0;
    label nTotal = 0;

    for (const label celli : cells_)
    {
        const scalar alphaI = alpha_[celli];

        nTotal++;

        if (alphaI < max_)
        {
            nLimit++;
        }
    }

    reduce(nLimit, sumOp<label>());
    reduce(nTotal, sumOp<label>());

    Info << nLimit << " cells are limited in total "<< nTotal << " cells." << endl;

    // handle boundaries in the case of 'all'
    if (selectionMode_ == smAll)
    {
        volScalarField::Boundary& Sbf = S.boundaryFieldRef();
        volScalarField::Boundary& alphaBf = alpha_.boundaryFieldRef();

        forAll(Sbf, patchi)
        {
            fvPatchScalarField& alphap = alphaBf[patchi];

            if (!alphap.fixesValue())
            {
                /*
                forAll(Sp, facei)
                {
                    if (Sp[facei] > maxScalar)
                    {
                        Sp[facei] = maxScalar;
                    }
                }
                */
            }
        }
    }

    // We've changed internal values so give
    // boundary conditions opportunity to correct
    S.correctBoundaryConditions();
}


// ************************************************************************* //
