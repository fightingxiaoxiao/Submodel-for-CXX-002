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

\*---------------------------------------------------------------------------*/

#include "gamma.H"
#include "addToRunTimeSelectionTable.H"
#include <random>
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace distributionModels
    {
        defineTypeNameAndDebug(gamma, 0);
        addToRunTimeSelectionTable(distributionModel, gamma, dictionary);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionModels::gamma::gamma
(
    const dictionary &dict, 
    Random& rndGen
)
    : distributionModel(typeName, dict, rndGen),
      alpha_(distributionModelDict_.get<scalar>("alpha")),
      beta_(distributionModelDict_.get<scalar>("beta")),
      meanValue_(distributionModelDict_.get<scalar>("meanValue")),
      maxValue_(distributionModelDict_.get<scalar>("maxValue")),
      minValue_(distributionModelDict_.get<scalar>("minValue"))
{
    check();
}

Foam::distributionModels::gamma::gamma(const gamma &p)
: 
    distributionModel(p),
    alpha_(p.alpha_),
    beta_(p.beta_),
    meanValue_(p.meanValue_),
    maxValue_(p.maxValue_),
    minValue_(p.minValue_)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributionModels::gamma::~gamma()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributionModels::gamma::sample() const
{
    std::random_device rd;
    std::mt19937 gen(rd());
    // C++标准库的gamma分布函数的形参虽然写作alpha,beta，但实际上是k,theta
    std::gamma_distribution<> gamma(alpha_, 1.0 / beta_);

    scalar rndSize = gamma(gen) / (alpha_ / beta_) * meanValue_;
    // 当尺寸超出上下界时重新生成随机尺寸
    while (rndSize > maxValue_ || rndSize < minValue_)
    {
        rndSize = gamma(gen) / (alpha_ / beta_) * meanValue_;
    }
    return rndSize;
}

// 计算颗粒团的最大、最小和平均粒径需要遍历整个颗粒云，这会增加不必要的计算开销，故此处返回值写得比较随意
// 需要仔细审视给出的参数是否合理，尤其是上下边界，需要尽可能地位于PDF边缘
Foam::scalar Foam::distributionModels::gamma::minValue() const
{
    return minValue_;
}

Foam::scalar Foam::distributionModels::gamma::maxValue() const
{
    return maxValue_;
}

Foam::scalar Foam::distributionModels::gamma::meanValue() const
{
    return meanValue_;
}

// ************************************************************************* //
