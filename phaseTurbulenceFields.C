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
    // Info << "The turbulence model is phase compressible :" << modelName_;

    // if (obr_.foundObject<compressible::turbulenceModel>(modelName_))
    // {
    //     return true;
    // }
    // else if (obr_.foundObject<incompressible::turbulenceModel>(modelName_))
    // {
    //     return false;
    // }
    // else if (obr_.foundObject<PhaseIncompressibleTurbulenceModel<transportModel>>(modelName_))
    // {
    //     Info << "The turbulence model is phase Incompressible :" << modelName_;
    //     return false;
    // }
    // else if (obr_.foundObject<PhaseCompressibleTurbulenceModel<transportModel>>(modelName_))
    // {
    //     Info << "The turbulence model is phase compressible :" << modelName_;
    //     return true;
    // }

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
            Info<< "storing fields:" << nl;
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


    const auto& model = obr_.lookupObject<turbulenceModel>(IOobject::groupName("turbulenceProperties", phaseName_));

    Info << "-------------------------------\n";
    // Info << "Model = " << model;
    Info << "Model2 = " << model;
    Info << "-------------------------------\n";
    
    // (
    //         turbulenceModel::propertiesName
    // );
    // Info << "The turbulence model is :" << model;
    // const transportModel& liquid = model->transport();
    // const twoPhaseSystem& fluid =
    //     refCast<const twoPhaseSystem>(liquid.fluid());
    // const transportModel& gas = fluid.otherPhase(liquid);

    // gasTurbulencePtr_ =
    // &U.db()
    // .lookupObject<PhaseCompressibleTurbulenceModel<transportModel>>
    // (
    //     IOobject::groupName
    //     (
    //         turbulenceModel::propertiesName,
    //         gas.name()
    //     )
    // );

    // const bool comp = compressible();

    // if (comp)
    // {
    //     // const auto& model =
    //     //     obr_.lookupObject<compressible::turbulenceModel>(modelName_);

        for (const word& f : fieldSet_)
        {
            switch (compressibleFieldNames_[f])
            {
                // case cfUPrime:
                // {
                //     // volVectorField uPrime(model.R().diag());
                //     processField<vector>(f, model.R().diag());
                //     break;
                // }
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
                // case cfLESRegion:
                // {
                //     auto* DESPtr = mesh_.cfindObject<DESModelBase>(modelName_);
                //     if (!DESPtr)
                //     {
                //         WarningInFunction
                //             << "Turbulence model is not a DES model - "
                //             << "skipping request for LESRegion" << endl;

                //         break;
                //     }

                    // processField<scalar>(f, DESPtr->LESRegion());
                    // break;
                // }
                // case cffd:
                // {
                //     auto* DESPtr = mesh_.cfindObject<DESModelBase>(modelName_);
                //     if (!DESPtr)
                //     {
                //         WarningInFunction
                //             << "Turbulence model is not a DES model - "
                //             << "skipping request for fd" << endl;

                //         break;
                //     }

                //     processField<scalar>(f, DESPtr->fd());
                //     break;
                // }
                default:
                {
                    FatalErrorInFunction
                        << "Invalid field selection" << abort(FatalError);
                }
            }
        }
    // }
    // else
    // {
    //     // const auto& model =
    //     //     obr_.lookupObject<incompressible::turbulenceModel>(modelName_);

    //     for (const word& f : fieldSet_)
    //     {
    //         switch (incompressibleFieldNames_[f])
    //         {
    //             case ifK:
    //             {
    //                 processField<scalar>(f, model.k());
    //                 break;
    //             }
    //             case ifEpsilon:
    //             {
    //                 processField<scalar>(f, model.epsilon());
    //                 break;
    //             }
    //             case ifOmega:
    //             {
    //                 processField<scalar>(f, model.omega());
    //                 break;
    //             }
    //             case ifNuTilda:
    //             {
    //                 processField<scalar>(f, nuTilda(model));
    //                 break;
    //             }
    //             case ifNut:
    //             {
    //                 processField<scalar>(f, model.nut());
    //                 break;
    //             }
    //             case ifNuEff:
    //             {
    //                 processField<scalar>(f, model.nuEff());
    //                 break;
    //             }
    //             case ifR:
    //             {
    //                 processField<symmTensor>(f, model.R());
    //                 break;
    //             }
    //             case ifDevReff:
    //             {
    //                 processField<symmTensor>(f, model.devReff());
    //                 break;
    //             }
    //             case ifL:
    //             {
    //                 processField<scalar>(f, L(model));
    //                 break;
    //             }
    //             case ifI:
    //             {
    //                 processField<scalar>(f, I(model));
    //                 break;
    //             }
    //             // case ifLESRegion:
    //             // {
    //             //     auto* DESPtr = mesh_.cfindObject<DESModelBase>(modelName_);
    //             //     if (!DESPtr)
    //             //     {
    //             //         WarningInFunction
    //             //             << "Turbulence model is not a DES model - "
    //             //             << "skipping request for LESRegion" << endl;

    //             //         break;
    //             //     }

    //             //     processField<scalar>(f, DESPtr->LESRegion());
    //             //     break;
    //             // }
    //             // case iffd:
    //             // {
    //             //     auto* DESPtr = mesh_.cfindObject<DESModelBase>(modelName_);
    //             //     if (!DESPtr)
    //             //     {
    //             //         WarningInFunction
    //             //             << "Turbulence model is not a DES model - "
    //             //             << "skipping request for fd" << endl;

    //             //         break;
    //             //     }

    //             //     processField<scalar>(f, DESPtr->fd());
    //             //     break;
    //             // }
    //             default:
    //             {
    //                 FatalErrorInFunction
    //                     << "Invalid field selection" << abort(FatalError);
    //             }
    //         }
    //     }
    // }

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
