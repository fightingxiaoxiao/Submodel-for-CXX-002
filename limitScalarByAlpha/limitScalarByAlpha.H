/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::fv::limitScalarByAlpha

Description
    Corrects scalar field (i.e. \c T) within a specified
    region by applying a given maximum scalar magnitude limit.

    Corrections applied to:
    \verbatim
      p/k/epslion/omega...    | Scalar
    \endverbatim

    Required fields:
    \verbatim
      p/k/epslion/omega...    | Scalar
    \endverbatim

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    limitScalarByAlpha1
    {
        // Mandatory entries (unmodifiable)
        type            limitScalarByAlpha;

        // Mandatory entries (runtime modifiable)
        max             <maxValue>;

        // Optional entries (runtime modifiable)
        scalarName               <Scalarname>;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                   | Type   | Reqd  | Dflt
      type      | Type name: limitScalarByAlpha      | word   | yes   | -
      max       | Maximum scalar limit [m/s]  | scalar | yes   | -
      Scalar    | Name of operand scalar field | word  | no    | Scalar
    \endtable

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink
     - \link cellSetOption.H \endlink

See also
  - Foam::fv::scalarDampingConstraint

SourceFiles
    limitScalarByAlpha.C

\*---------------------------------------------------------------------------*/

#ifndef limitScalarByAlpha_H
#define limitScalarByAlpha_H

#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class limitScalarByAlpha Declaration
\*---------------------------------------------------------------------------*/

class limitScalarByAlpha
:
    public cellSetOption
{
protected:

    // Protected Data

        //- Name of operand scalar field
        word scalarName_;
        
        //- Maximum scalar magnitude
        scalar max_;

        word phaseName_;


public:

    //- Runtime type information
    TypeName("limitScalarByAlpha");


    // Constructors

        //- Construct from components
        limitScalarByAlpha
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        limitScalarByAlpha(const limitScalarByAlpha&) = delete;

        //- No copy assignment
        void operator=(const limitScalarByAlpha&) = delete;


    //- Destructor
    virtual ~limitScalarByAlpha() = default;


    // Member Functions

        //- Read dictionary
        virtual bool read(const dictionary& dict);

        //- Correct the scalar field
        virtual void correct(volScalarField& scalarField);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
