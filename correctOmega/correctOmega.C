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

#include "correctOmega.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(correctOmega, 0);
    addToRunTimeSelectionTable(option, correctOmega, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::correctOmega::correctOmega
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    omegaName_(coeffs_.get<word>("omegaName")),
    alphaName_(coeffs_.get<word>("alphaName")),
    nutName_(coeffs_.get<word>("nutName")),
    lowerLimit_(coeffs_.get<scalar>("lowerLimit")),
    upperLimit_(coeffs_.get<scalar>("upperLimit")),
    B_(coeffs_.get<scalar>("B")),
    betaStar_(coeffs_.getOrDefault<scalar>("betaStar", 0.09))
{
    fieldNames_.setSize(1, omegaName_);
    applied_.setSize(1, false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::correctOmega::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.readEntry("omegaName", omegaName_);
        coeffs_.readEntry("alphaName", alphaName_);
        coeffs_.readEntry("nutName", nutName_);
        coeffs_.readEntry("lowerLimit", lowerLimit_);
        coeffs_.readEntry("upperLimit", upperLimit_);
        coeffs_.readEntry("B", B_);
        coeffs_.readEntry("betaStar", betaStar_);
        return true;
    }

    return false;
}


void Foam::fv::correctOmega::correct(volScalarField& S)
{
    const volScalarField &alpha_ = mesh_.lookupObject<volScalarField>(alphaName_);
    const volScalarField &nut_ = mesh_.lookupObject<volScalarField>(nutName_);

    label nCorrect = 0;
    label nTotal = 0;

    for (const label celli : cells_)
    {
        const scalar alphaI = alpha_[celli];

        const scalar deltaX = cbrt(mesh_.V()[celli]); 
        
        if (alphaI < upperLimit_ && alphaI > lowerLimit_)
        {
            S[celli] = B_*6.*nut_[celli]/betaStar_/deltaX/deltaX;
            nCorrect++;
        }

        nTotal++;
    }
    
    reduce(nCorrect, sumOp<label>());
    reduce(nTotal, sumOp<label>());

    Info <<"multiphaseOmegaCorrect: "<< nCorrect << " cells are corrected in total "<< nTotal << " cells." << endl;

    // handle boundaries in the case of 'all'
    if (selectionMode_ == smAll)
    {
        //volScalarField::Boundary& Sbf = S.boundaryFieldRef();
        //const volScalarField::Boundary& alphaBf = alpha_.boundaryField();
        /*
        forAll(Sbf, patchi)
        {
            //const fvPatchScalarField& alphap = alphaBf[patchi];
            //fvPatchScalarField& Sp = Sbf[patchi];

            if (!Sp.fixesValue())
            {
                
                forAll(Sp, facei)
                {
                    if (alphap[facei] < maxAlpha)
                    {
                        Sp[facei] = 0.;
                        Sp[facei] = B_*6.*nut_[celli]/betaStar_/deltaX/delatX;
                    }
                }
                
            }
        }
        */
    }

    // We've changed internal values so give
    // boundary conditions opportunity to correct
    S.correctBoundaryConditions();
}


// ************************************************************************* //
