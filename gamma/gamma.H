/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::distributionModels::gamma

Description
    gamma distribution model

SourceFiles
    gamma.C

\*---------------------------------------------------------------------------*/

#ifndef distributionModels_gamma_H
#define distributionModels_gamma_H

#include "distributionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace distributionModels
{

/*---------------------------------------------------------------------------*\
                        Class gamma Declaration
\*---------------------------------------------------------------------------*/

class gamma
:
    public distributionModel
{
    // Private data

        scalar alpha_;

        scalar beta_;

        // 平均尺寸
        scalar meanValue_;

        // 用于截断最大尺寸
        scalar maxValue_;

        // 用于截断最小尺寸
        scalar minValue_;

public:

    //- Runtime type information
    TypeName("gamma");


    // Constructors

        //- Construct from components
        gamma(const dictionary& dict, Random& rndGen);

        //- Construct copy
        gamma(const gamma& p);

        //- Construct and return a clone
        virtual autoPtr<distributionModel> clone() const
        {
            return autoPtr<distributionModel>(new gamma(*this));
        }


    //- Destructor
    virtual ~gamma();


    // Member Functions

        //- Sample the distributionModel
        virtual scalar sample() const;

        //- Return the minimum value
        virtual scalar minValue() const;

        //- Return the maximum value
        virtual scalar maxValue() const;

        //- Return the mean value
        virtual scalar meanValue() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace distributionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
