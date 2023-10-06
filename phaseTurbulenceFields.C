/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "phaseTurbulenceFields.H"
#include "turbulentTransportModel.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "turbulentFluidThermoModel.H"
#include "DESModelBase.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(phaseTurbulenceFields, 0);
    addToRunTimeSelectionTable(functionObject, phaseTurbulenceFields, dictionary);
}
}

const Foam::Enum
<
    Foam::functionObjects::phaseTurbulenceFields::compressibleField
>
Foam::functionObjects::phaseTurbulenceFields::compressibleFieldNames_
({
    { compressibleField::cfK, "k" },
    { compressibleField::cfEpsilon, "epsilon" },
    { compressibleField::cfOmega, "omega" },
    { compressibleField::cfNuTilda, "nuTilda" },
    { compressibleField::cfMut, "mut" },
    { compressibleField::cfMuEff, "muEff" },
    { compressibleField::cfAlphat, "alphat" },
    { compressibleField::cfAlphaEff, "alphaEff" },
    { compressibleField::cfR, "R" },
    { compressibleField::cfDevRhoReff, "devRhoReff" },
    { compressibleField::cfL, "L" },
    { compressibleField::cfI, "I" },
    { compressibleField::cfLESRegion, "LESRegion" },
    { compressibleField::cffd, "fd" },
});

//TODO : remove these
const Foam::Enum
<
    Foam::functionObjects::phaseTurbulenceFields::incompressibleField
>
Foam::functionObjects::phaseTurbulenceFields::incompressibleFieldNames_
({
    { incompressibleField::ifK, "k" },
    { incompressibleField::ifEpsilon, "epsilon" },
    { incompressibleField::ifOmega, "omega" },
    { incompressibleField::ifNuTilda, "nuTilda" },
    { incompressibleField::ifNut, "nut" },
    { incompressibleField::ifNuEff, "nuEff" },
    { incompressibleField::ifR, "R" },
    { incompressibleField::ifDevReff, "devReff" },
    { incompressibleField::ifL, "L" },
    { incompressibleField::ifI, "I" },
    { incompressibleField::ifLESRegion, "LESRegion" },
    { incompressibleField::iffd, "fd" },
});



// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::phaseTurbulenceFields::initialise()
{


    for (const word& f : fieldSet_)
    {
        const word localName(IOobject::scopedName(prefix_,  IOobject::groupName(f,phaseName_)));
        // const word localName(IOobject::scopedName(prefix_, f));

        if (obr_.found(localName))
        {
            WarningInFunction
                << "Cannot store turbulence field " << localName
                << " since an object with that name already exists"
                << nl << endl;

            fieldSet_.unset(f);
        }
    }

    initialised_ = true;
}


bool Foam::functionObjects::phaseTurbulenceFields::compressible()
{
    //TODO : remove this no longer used function
    FatalErrorInFunction
        << "Turbulence model not found in database, deactivating"
        << exit(FatalError);

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::phaseTurbulenceFields::phaseTurbulenceFields
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    initialised_(false),
    prefix_(dict.getOrDefault<word>("prefix", "turbulenceProperties")),
    phaseName_(dict.get<word>("phase")),
    fieldSet_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::phaseTurbulenceFields::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        dict.readIfPresent("prefix", prefix_);

        if (dict.found("field"))
        {
            fieldSet_.insert(dict.get<word>("field"));
        }
        else
        {
            fieldSet_.insert(dict.get<wordList>("fields"));
        }

        Info<< type() << " " << name() << ": ";
        if (fieldSet_.size())
        {
            Info<< "Storing fields:" << nl;
            for (const word& f : fieldSet_)
            {
                Info<< "    " << IOobject::scopedName(prefix_,  IOobject::groupName(f,phaseName_)) << nl;
            }
            Info<< endl;
        }
        else
        {
            Info<< "no fields requested to be stored" << nl << endl;
        }

        initialised_ = false;

        return true;
    }

    return false;
}


bool Foam::functionObjects::phaseTurbulenceFields::execute()
{
    if (!initialised_)
    {
        initialise();
    }

    // get the current phase turbulence model used
    // TODO : rename the enums to reflect the case that we dont segregate based on the model used
    //! multiphase models are virtually all compressible models
    const auto& model = obr_.lookupObject<turbulenceModel>(IOobject::groupName("turbulenceProperties", phaseName_));

        for (const word& f : fieldSet_)
        {
            switch (compressibleFieldNames_[f])
            {
                // TODO: add a way to calculate the rsm fluctuations directly UPrime (or the variance UPrime2)
                // const volScalarField uPrime(sqrt((2.0/3.0)*model.k()));
                case cfK:
                {
                    processField<scalar>(f, model.k());
                    break;
                }
                case cfEpsilon:
                {
                    processField<scalar>(f, model.epsilon());
                    break;
                }
                case cfOmega:
                {
                    processField<scalar>(f, model.omega());
                    break;
                }
                case cfNuTilda:
                {
                    processField<scalar>(f, nuTilda(model));
                    break;
                }
                case cfMut:
                {
                    processField<scalar>(f, model.mut());
                    break;
                }
                case cfMuEff:
                {
                    processField<scalar>(f, model.muEff());
                    break;
                }
                case cfAlphat:
                {
                    // processField<scalar>(f, model.alphat());
                    break;
                }
                case cfAlphaEff:
                {
                    // processField<scalar>(f, model.alphaEff());
                    break;
                }
                case cfR:
                {
                    processField<symmTensor>(f, model.R());
                    break;
                }
                case cfDevRhoReff:
                {
                    // processField<symmTensor>(f, model.devRhoReff());
                    break;
                }
                case cfL:
                {
                    processField<scalar>(f, L(model));
                    break;
                }
                case cfI:
                {
                    processField<scalar>(f, I(model));
                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Invalid field selection" << abort(FatalError);
                }
            }
    return true;
}


bool Foam::functionObjects::phaseTurbulenceFields::write()
{
    for (const word& f : fieldSet_)
    {
        // const word localName(IOobject::scopedName(prefix_, f));
        const word localName(IOobject::scopedName(prefix_, IOobject::groupName(f,phaseName_)));
        Info << "Writing " << localName;
        writeObject(localName);
        Info << "\n";
    }
    Info<< endl;

    return true;
}


// ************************************************************************* //
