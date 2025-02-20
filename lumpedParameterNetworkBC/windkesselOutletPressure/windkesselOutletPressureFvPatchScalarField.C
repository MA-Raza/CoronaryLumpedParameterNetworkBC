/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "windkesselOutletPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "backwardDdtScheme.H"
#include "word.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// (No additional details required for static members here)

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Constructor: Initialize with default values
Foam::windkesselOutletPressureFvPatchScalarField::
windkesselOutletPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF), // Call base class constructor
    Rp_(1),     // Proximal resistance
    Rd_(1),     // Distal resistance
    C_(1),      // Compliance
    L_(1),      // Inductance
    Pd_(0),     // Distal pressure
    rho_(1),    // Fluid density
    Pooo_(0), Poo_(0), Po_(0), Pn_(0), // Initialize pressures to zero
    Qooo_(0), Qoo_(0), Qo_(0), Qn_(0), // Initialize flow rates to zero
    windkesselModel_(WK3),     // Default windkessel model
    diffScheme_(secondOrder),  // Default differencing scheme
    timeIndex_(-1)  // Initialize time index to invalid
{}

//- Constructor: Initialize from dictionary (typically used in user input)
Foam::windkesselOutletPressureFvPatchScalarField::
windkesselOutletPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),      // Call base class constructor
    Rp_(dict.lookupOrDefault<scalar>("Rp", 1)),     // Proximal resistance
    Rd_(dict.lookupOrDefault<scalar>("Rd", 1)),     // Distal resistance
    C_(dict.lookupOrDefault<scalar>("C", 1)),       // Compliance
    L_(dict.lookupOrDefault<scalar>("L", 1)),       // Inductance
    Pd_(dict.lookupOrDefault<scalar>("Pd", 0)),     // Distal pressure
    rho_(dict.lookupOrDefault<scalar>("rho", 1)),   // Fluid density
    Pooo_(dict.lookupOrDefault<scalar>("Pooo", 0)),
    Poo_(dict.lookupOrDefault<scalar>("Poo", 0)),
    Po_(dict.lookupOrDefault<scalar>("Po", 0)),
    Pn_(dict.lookupOrDefault<scalar>("Pn", 0)),
    Qooo_(dict.lookupOrDefault<scalar>("Qooo", 0)),
    Qoo_(dict.lookupOrDefault<scalar>("Qoo", 0)),
    Qo_(dict.lookupOrDefault<scalar>("Qo", 0)),
    Qn_(dict.lookupOrDefault<scalar>("Qn", 0)),
    timeIndex_(-1)  // Time index reset
{

    Info << "\n\nApplying windkesselOutletPressure BC on patch: " << patch().name() << endl;

    //- Retrieve the model as a string and map to the enumerator
    word WKModelStr = dict.lookupOrDefault<word>("windkesselModel", "WK3");

    if (WKModelStr == "Resistive")
    {
        windkesselModel_ = Resistive;

        //- Print Windkessel model properties for debugging/verification
        Info << "\n\nProperties of Windkessel Model: \n"
             << "Model Type: 1-Element (Resistive) \n"
             << "Distal Resistance: Rd = " << Rd_ << " kgm^-4s^-1 \n"
             << "Distal Pressure: Pd = " << Pd_ << " Pa \n"
             << "Density: rho = " << rho_ << " kgm^-3 \n" << endl;
    }
    else if (WKModelStr == "WK2")
    {
        windkesselModel_ = WK2;

        //- Print Windkessel model properties for debugging/verification
        Info << "\n\nProperties of Windkessel Model: \n"
             << "Model Type: 2-Element (RC) \n"
             << "Distal Resistance: Rd = " << Rd_ << " kgm^-4s^-1 \n"
             << "Compliance: C = " << Rd_ << " m^4s^2kg^-1 \n"
             << "Distal Pressure: Pd = " << C_ << " Pa \n"
             << "Density: rho = " << rho_ << " kgm^-3 \n" << endl;
    }
    else if (WKModelStr == "WK3")
    {
        windkesselModel_ = WK3;

        //- Print Windkessel model properties for debugging/verification
        Info << "\n\nProperties of Windkessel Model: \n"
             << "Model Type: 3-Element (RCR) \n"
             << "Proximal Resistance: Rp = " << Rp_ << " kgm^-4s^-1 \n"
             << "Distal Resistance: Rd = " << Rd_ << " kgm^-4s^-1 \n"
             << "Compliance: C = " << C_ << " m^4s^2kg^-1 \n"
             << "Distal Pressure: Pd = " << Pd_ << " Pa \n"
             << "Density: rho = " << rho_ << " kgm^-3 \n" << endl;
    }
    else if (WKModelStr == "WK4Series")
    {
        windkesselModel_ = WK4Series;
        
        //- Print Windkessel model properties for debugging/verification
        Info << "\n\nProperties of Windkessel Model: \n"
             << "Model Type: 4-Element (RCRL with L in Series to Rp) \n"
             << "Proximal Resistance: Rp = " << Rp_ << " kgm^-4s^-1 \n"
             << "Distal Resistance: Rd = " << Rd_ << " kgm^-4s^-1 \n"
             << "Compliance: C = " << C_ << " m^4s^2kg^-1 \n"
             << "Inertance: L = " << L_ << " kgm^-4 \n"
             << "Distal Pressure: Pd = " << Pd_ << " Pa \n"
             << "Density: rho = " << rho_ << " kgm^-3 \n" << endl;
    }
    else if (WKModelStr == "WK4Parallel")
    {
        windkesselModel_ = WK4Parallel;

        //- Print Windkessel model properties for debugging/verification
        Info << "\n\nProperties of Windkessel Model: \n"
             << "Model Type: 4-Element (RCRL with L in Parallel to Rp) \n"
             << "Proximal Resistance: Rp = " << Rp_ << " kgm^-4s^-1 \n"
             << "Distal Resistance: Rd = " << Rd_ << " kgm^-4s^-1 \n"
             << "Compliance: C = " << C_ << " m^4s^2kg^-1 \n"
             << "Inertance: L = " << L_ << " kgm^-4 \n"
             << "Distal Pressure: Pd = " << Pd_ << " Pa \n"
             << "Density: rho = " << rho_ << " kgm^-3 \n" << endl;
    }
    else
    {
        FatalErrorInFunction << "\n\nUnknown Windkessel Model: " << WKModelStr
                             << "\nValid Windkessel Models (windkesselModel) are : \n\n"
                             << " 5 \n ( \n Resistive \n WK2 \n WK3 \n WK4Series \n WK4Parallel \n ) \n"
                             << exit(FatalError);
    }

    //- Retrieve the scheme as a string and map to the enumerator
    word schemeStr = dict.lookupOrDefault<word>("diffScheme", "secondOrder");

    Info << "Differencing Scheme for Windkessel Model: " << schemeStr << endl;

    if (schemeStr == "firstOrder")
    {
        diffScheme_ = firstOrder;
    }
    else if (schemeStr == "secondOrder")
    {
        diffScheme_ = secondOrder;
    }
    else
    {
        FatalErrorInFunction << "\n\nUnknown Differencing Scheme: " << schemeStr
                             << "\nValid Differencing Schemes (diffScheme) are : \n\n"
                             << " 2 \n ( \n firstOrder \n secondOrder \n ) \n"
                             << exit(FatalError);
    }

}

//- Constructor: Map an existing field onto a new patch
Foam::windkesselOutletPressureFvPatchScalarField::
windkesselOutletPressureFvPatchScalarField
(
    const windkesselOutletPressureFvPatchScalarField& ptf, // Existing field
    const fvPatch& p, // New patch to map onto
    const DimensionedField<scalar, volMesh>& iF, // Internal field reference
    const fvPatchFieldMapper& mapper // Field mapper for mapping operations
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper), // Map base class properties
    Rp_(ptf.Rp_),     // Copy proximal resistance parameter
    Rd_(ptf.Rd_),     // Copy distal resistance parameter
    C_(ptf.C_),       // Copy compliance parameter
    L_(ptf.L_),       // Copy inertance parameter
    Pd_(ptf.Pd_),     // Copy distal pressure parameter
    rho_(ptf.rho_),   // Copy fluid density
    Pooo_(ptf.Pooo_), // Copy previous state pressure variables
    Poo_(ptf.Poo_),
    Po_(ptf.Po_),
    Pn_(ptf.Pn_),
    Qooo_(ptf.Qooo_), // Copy previous state flow rate variables
    Qoo_(ptf.Qoo_),
    Qo_(ptf.Qo_),
    Qn_(ptf.Qn_),
    windkesselModel_(ptf.windkesselModel_), // Copy Windkessel model type
    diffScheme_(ptf.diffScheme_), // Copy numerical differencing scheme
    timeIndex_(ptf.timeIndex_) // Copy time index
{}

//- Copy constructor: Create a duplicate field with all properties
Foam::windkesselOutletPressureFvPatchScalarField::
windkesselOutletPressureFvPatchScalarField
(
    const windkesselOutletPressureFvPatchScalarField& wkpsf // Source field
)
:
    fixedValueFvPatchScalarField(wkpsf), // Copy base class properties
    Rp_(wkpsf.Rp_),     // Copy proximal resistance
    Rd_(wkpsf.Rd_),     // Copy distal resistance
    C_(wkpsf.C_),       // Copy compliance
    L_(wkpsf.L_),       // Copy inertance
    Pd_(wkpsf.Pd_),     // Copy distal pressure
    rho_(wkpsf.rho_),   // Copy fluid density
    Pooo_(wkpsf.Pooo_), // Copy pressure variables
    Poo_(wkpsf.Poo_),
    Po_(wkpsf.Po_),
    Pn_(wkpsf.Pn_),
    Qooo_(wkpsf.Qooo_), // Copy flow rate variables
    Qoo_(wkpsf.Qoo_),
    Qo_(wkpsf.Qo_),
    Qn_(wkpsf.Qn_),
    windkesselModel_(wkpsf.windkesselModel_), // Copy Windkessel model type
    diffScheme_(wkpsf.diffScheme_), // Copy differencing scheme
    timeIndex_(wkpsf.timeIndex_) // Copy time index    
{}

//- Copy constructor with new internal field reference
Foam::windkesselOutletPressureFvPatchScalarField::
windkesselOutletPressureFvPatchScalarField
(
    const windkesselOutletPressureFvPatchScalarField& wkpsf, // Source field
    const DimensionedField<scalar, volMesh>& iF // New internal field
)
:
    fixedValueFvPatchScalarField(wkpsf, iF), // Map base class properties
    Rp_(wkpsf.Rp_),     // Copy proximal resistance
    Rd_(wkpsf.Rd_),     // Copy distal resistance
    C_(wkpsf.C_),       // Copy compliance
    L_(wkpsf.L_),       // Copy inertance
    Pd_(wkpsf.Pd_),     // Copy distal pressure
    rho_(wkpsf.rho_),   // Copy fluid density
    Pooo_(wkpsf.Pooo_), // Copy pressure variables
    Poo_(wkpsf.Poo_),
    Po_(wkpsf.Po_),
    Pn_(wkpsf.Pn_),
    Qooo_(wkpsf.Qooo_), // Copy flow rate variables
    Qoo_(wkpsf.Qoo_),
    Qo_(wkpsf.Qo_),
    Qn_(wkpsf.Qn_),
    windkesselModel_(wkpsf.windkesselModel_), // Copy Windkessel model type
    diffScheme_(wkpsf.diffScheme_), // Copy differencing scheme
    timeIndex_(wkpsf.timeIndex_) // Copy time index
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Update coefficients for boundary condition
void Foam::windkesselOutletPressureFvPatchScalarField::updateCoeffs()
{
    //- Check if coefficients have already been updated for this time step
    if (updated())
    {
        return;
    }

    //- Retrieve the time step size from the database
    const scalar dt = db().time().deltaTValue();

    //- Retrieve the flux field (phi) from the database
    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>("phi");

    //- Retrieve the boundary pressure field from the database
    const fvPatchField<scalar>& p = 
        db().lookupObject<volScalarField>("p").boundaryField()[this->patch().index()];

    //- Calculate the total area of the patch
    scalar area = gSum(patch().magSf());

    //- Update pressure and flow rate histories if the time index has advanced
    if (db().time().timeIndex() != timeIndex_)
    {
        timeIndex_ = db().time().timeIndex();

        //- Update pressure history variables
        Pooo_ = Poo_;
        Poo_ = Po_;
        Po_ = Pn_;
        Pn_ = gSum(p * patch().magSf()) / area; // Compute mean pressure over the patch

        //- Update flow rate history variables
        Qooo_ = Qoo_;
        Qoo_ = Qo_;
        Qo_ = Qn_;
        Qn_ = gSum(phi.boundaryField()[patch().index()]); // Compute total flux
    }

    //- Calculate new pressure using backward differencing scheme
    if
    (
        db().time().timeIndex()>1 // Perform calculations only if the time index is greater than 1
    )
    {
        //- Switch to the specified windkessel model
        switch (windkesselModel_)
        {
            case Resistive: // Single-element resistive model (R)

                Pn_ = Rd_ * Qn_ + Pd_;

                break;
            
            case WK2: // 2-element windkessel model (RC)

                //- Switch to the specified differencing scheme
                switch (diffScheme_)
                {
                    case firstOrder:

                        Pn_ = Rd_ * Qn_
                            + Pd_
                            + Rd_ * C_ * Po_ / dt;

                        Pn_ /= (1.0 + Rd_ * C_ / dt) + SMALL;

                        break;

                    case secondOrder:

                        Pn_ = Rd_ * Qn_
                            + Pd_
                            - Rd_ * C_ * ((Poo_ - 4 * Po_) / (2 * dt));

                        Pn_ /= (1.0 + 3 * Rd_ * C_ / (2 * dt)) + SMALL;

                        break;

                    default:
                    
                        FatalErrorInFunction << "Unknown differencing scheme!" << exit(FatalError);
                }

                break;

            case WK3: // 3-element windkessel model (RCR)

                //- Switch to the specified differencing scheme
                switch (diffScheme_)
                {
                    case firstOrder:

                        Pn_ = Rp_ * Rd_ * C_ * ((Qn_ - Qo_) / dt)
                            + (Rp_ + Rd_) * Qn_
                            + Pd_
                            + Rd_ * C_ * (Po_ / dt);

                        Pn_ /= (1.0 + Rd_ * C_ / dt) + SMALL;

                        break;

                    case secondOrder:

                        Pn_ = Rp_ * Rd_ * C_ * ((3 * Qn_ - 4 * Qo_ + Qoo_) / (2 * dt))
                            + (Rp_ + Rd_) * Qn_
                            + Pd_
                            - Rd_ * C_ * ((Poo_ - 4 * Po_) / (2 * dt));

                        Pn_ /= (1.0 + 3 * Rd_ * C_  / (2 * dt)) + SMALL;

                        break;

                    default:
                    
                        FatalErrorInFunction << "Unknown differencing scheme!" << exit(FatalError);
                }

                break;
            
            case WK4Series: // 4-element series windkessel model (RCRL-Series)

                //- Switch to the specified differencing scheme
                switch (diffScheme_)
                {
                    case firstOrder:

                        Pn_ = (Rp_ + Rd_) * Qn_
                            + (L_ + Rp_ * Rd_ * C_) * ((Qn_ - Qo_) / dt)
                            + Rd_ * C_ * L_ * ((Qn_ - 2 * Qo_ + Qoo_) / (pow (dt, 2)))
                            + Pd_
                            + Rd_ * C_ * (Po_ / dt);

                        Pn_ /= (1.0 + Rd_ * C_ / dt) + SMALL;

                        break;

                    case secondOrder:

                        Pn_ = (Rp_ + Rd_) * Qn_
                            + (L_ + Rp_ * Rd_ * C_) * ((3 * Qn_ - 4 * Qo_ + Qoo_) / (2 * dt))
                            + Rd_ * C_ * L_ * ((2 * Qn_ - 5 * Qo_ + 4 * Qoo_ - Qooo_) / (pow (dt, 2)))
                            + Pd_
                            - Rd_ * C_ * ((Poo_ - 4 * Po_) / (2 * dt));

                        Pn_ /= (1.0 + 3 * Rd_ * C_  / (2 * dt)) + SMALL;

                        break;

                    default:
                    
                        FatalErrorInFunction << "Unknown differencing scheme!" << exit(FatalError);
                }

                break;

            case WK4Parallel: // 4-element parallel windkessel model (RCRL-Parallel)

                //- Switch to the specified differencing scheme
                switch (diffScheme_)
                {
                    case firstOrder:

                        Pn_ = Rd_ * Qn_
                            + L_ * (1 + Rd_ / Rp_ ) * ((Qn_ - Qo_) / dt)
                            + Rd_ * C_ * L_ * ((Qn_ - 2 * Qo_ + Qoo_) / (pow (dt, 2)))
                            + Pd_
                            + ((L_ + Rp_ * Rd_ * C_) / Rp_) * (Po_ / dt)
                            + (Rd_ * C_ * L_ / Rp_) * ((2 * Po_ - Poo_) / (pow (dt, 2)));

                        Pn_ /= (1.0 + ((L_ + Rp_ * Rd_ * C_) / (Rp_ * dt)) + (Rd_ * C_ * L_ / (Rp_ * pow (dt, 2)))) + SMALL;

                        break;

                    case secondOrder:

                        Pn_ = Rd_ * Qn_
                            + L_ * (1 + Rd_ / Rp_ ) * ((3 * Qn_ - 4 * Qo_ + Qoo_) / (2 * dt))
                            + Rd_ * C_ * L_ * ((2 * Qn_ - 5 * Qo_ + 4 * Qoo_ - Qooo_) / (pow (dt, 2)))
                            + Pd_
                            - ((L_ + Rp_ * Rd_ * C_) / Rp_) * ((Poo_ - 4 * Po_) / (2 * dt))
                            - (Rd_ * C_ * L_ / Rp_) * ((-5 * Po_ + 4 * Poo_ - Pooo_) / (pow (dt, 2)));

                        Pn_ /= (1.0 + (3 * (L_ + Rp_ * Rd_ * C_) / (2 * Rp_ * dt)) + (2 * Rd_ * C_ * L_ / (Rp_ * pow (dt, 2)))) + SMALL;

                        break;

                    default:
                    
                        FatalErrorInFunction << "Unknown differencing scheme!" << exit(FatalError);
                }

                break;

            default:

                FatalErrorInFunction << "Unknown Windkessel Model!" << exit(FatalError);
        }

    }

    //- Apply implicit pressure update to the boundary field
    operator==(Pn_ / (rho_ + SMALL));

    //- Call base class function to finalize the update
    fixedValueFvPatchScalarField::updateCoeffs();
}

//- Write the boundary condition properties to an output stream
void Foam::windkesselOutletPressureFvPatchScalarField::write(Ostream& os) const
{
    //- Call the base class function to write common boundary field properties
    fvPatchScalarField::write(os);

    //- Write the common windkessel model details to the output stream
    os.writeKeyword("Rd") << Rd_ << token::END_STATEMENT << nl; // Resistance parameter
    os.writeKeyword("Pd") << Pd_ << token::END_STATEMENT << nl; // Prescribed pressure
    os.writeKeyword("rho") << rho_ << token::END_STATEMENT << nl; // Fluid density

    //- Write the specific windkessel model details to the output stream
    switch (windkesselModel_)
    {
        case Resistive:
            os.writeKeyword("windkesselModel") << "Resistive" << token::END_STATEMENT << nl;
            break;

        case WK2:
            os.writeKeyword("windkesselModel") << "WK2" << token::END_STATEMENT << nl;
            os.writeKeyword("C") << C_ << token::END_STATEMENT << nl; // Compliance
            break;

        case WK3:
            os.writeKeyword("windkesselModel") << "WK3" << token::END_STATEMENT << nl;
            os.writeKeyword("Rp") << Rp_ << token::END_STATEMENT << nl; // Proximal resistance
            os.writeKeyword("C") << C_ << token::END_STATEMENT << nl; // Compliance
            break;

        case WK4Series:
            os.writeKeyword("windkesselModel") << "WK4Series" << token::END_STATEMENT << nl;
            os.writeKeyword("Rp") << Rp_ << token::END_STATEMENT << nl;
            os.writeKeyword("C") << C_ << token::END_STATEMENT << nl;
            os.writeKeyword("L") << L_ << token::END_STATEMENT << nl; // Inductance
            break;

        case WK4Parallel:
            os.writeKeyword("windkesselModel") << "WK4Parallel" << token::END_STATEMENT << nl;
            os.writeKeyword("Rp") << Rp_ << token::END_STATEMENT << nl;
            os.writeKeyword("C") << C_ << token::END_STATEMENT << nl;
            os.writeKeyword("L") << L_ << token::END_STATEMENT << nl;
            break;

        default:
            FatalErrorInFunction << "Unknown Windkessel Model: " << windkesselModel_ << exit(FatalError);
    }

    //- Write the differencing scheme used for calculations
    switch (diffScheme_)
    {
        case firstOrder:
            os.writeKeyword("diffScheme") << "firstOrder" << token::END_STATEMENT << nl;
            break;

        case secondOrder:
            os.writeKeyword("diffScheme") << "secondOrder" << token::END_STATEMENT << nl;
            break;

        default:
            FatalErrorInFunction << "Unknown differencing scheme: " << diffScheme_ << exit(FatalError);
    }

    //- Write flow history parameters
    os.writeKeyword("Qooo") << Qooo_ << token::END_STATEMENT << nl; // Third past flow rate
    os.writeKeyword("Qoo") << Qoo_ << token::END_STATEMENT << nl; // Second past flow rate
    os.writeKeyword("Qo") << Qo_ << token::END_STATEMENT << nl; // Previous flow rate
    os.writeKeyword("Qn") << Qn_ << token::END_STATEMENT << nl; // Current flow rate

    //- Write pressure history parameters
    os.writeKeyword("Pooo") << Pooo_ << token::END_STATEMENT << nl; // Third past pressure
    os.writeKeyword("Poo") << Poo_ << token::END_STATEMENT << nl; // Second past pressure
    os.writeKeyword("Po") << Po_ << token::END_STATEMENT << nl; // Previous pressure
    os.writeKeyword("Pn") << Pn_ << token::END_STATEMENT << nl; // Current pressure

    //- Write the boundary field value entry to the output stream
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    //- Define the boundary field type for this class
    makePatchTypeField
    (
        fvPatchScalarField,                         // Base class
        windkesselOutletPressureFvPatchScalarField  // Derived class
    );
}

// ************************************************************************* //
