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


// Include the header file defining the class and related dependencies
#include "coronaryOutletPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H" // Macro for runtime selection table
#include "fvPatchFieldMapper.H"         // Field mapper utility
#include "volFields.H"                  // Volume fields
#include "surfaceFields.H"              // Surface fields
#include "scalar.H"                     // Scalar type definition
#include <fstream>                      // File handling for reading/writing
#include <string.h>                     // String manipulation utilities
#include "IOdictionary.H"               // Dictionary I/O operations
#include "fileName.H"                   // FileName handling
#include "OSspecific.H"                 // OS-specific operations (e.g., paths)
#include <cstdlib>                      // Environment variable handling


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// (None in this section)

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Default constructor: Initializes the class with default values for all parameters.
Foam::coronaryOutletPressureFvPatchScalarField::
coronaryOutletPressureFvPatchScalarField
(
    const fvPatch& p,                                 // Patch reference
    const DimensionedField<Foam::scalar, volMesh>& iF // Internal field reference
)
:
    fixedValueFvPatchScalarField(p, iF), // Call the base class constructor
    Ra_(1),                              // Initialize default arterial resistance
    Ram_(1),                             // Initialize micro-arterial resistance
    Rv_(1),                              // Initialize venous resistance
    Rvm_(1),                             // Initialize micro-venous resistance
    Ca_(1),                              // Initialize arterial compliance
    Cim_(1),                             // Initialize intramyocardial compliance
    PimScaling_(1),                      // Initialize scaling factor for intramyocardial pressure
    Pv_(0),                              // Default distal pressure
    rho_(1),                             // Default fluid density
    Qoo_(0), Qo_(0), Qn_(0),             // Initialize flow rate history
    Poo_(0), Po_(0), Pn_(0),             // Initialize pressure history
    Pioo_(0), Pio_(0), Pin_(0),          // Initialize intermediate pressure history
    Pimoo_(0), Pimo_(0), Pimn_(0),       // Initialize intramyocardial pressure history
    diffScheme_(secondOrder),            // Use second-order differencing by default
    PimFile_("PimData"),                   //Default name of the Pim data file
    timeIndex_(-1)                       // Set an invalid time index as the initial state
{}

//- Constructor from a dictionary, typically used for user input configuration.
Foam::coronaryOutletPressureFvPatchScalarField::
coronaryOutletPressureFvPatchScalarField
(
    const fvPatch& p,                                 // Patch reference
    const DimensionedField<Foam::scalar, volMesh>& iF, // Internal field reference
    const dictionary& dict                           // Dictionary containing user-provided values
)
:
    fixedValueFvPatchScalarField(p, iF, dict), // Call the base class constructor
    Ra_(dict.lookupOrDefault<Foam::scalar>("Ra", 1)),       // Arterial resistance
    Ram_(dict.lookupOrDefault<Foam::scalar>("Ram", 1)),     // Micro-arterial resistance
    Rv_(dict.lookupOrDefault<Foam::scalar>("Rv", 1)),       // Venous resistance
    Rvm_(dict.lookupOrDefault<Foam::scalar>("Rvm", 1)),     // Micro-venous resistance
    Ca_(dict.lookupOrDefault<Foam::scalar>("Ca", 1)),       // Arterial compliance
    Cim_(dict.lookupOrDefault<Foam::scalar>("Cim", 1)),     // Intramyocardial compliance
    PimScaling_(dict.lookupOrDefault<Foam::scalar>("PimScaling", 1)), // Pim scaling factor
    Pv_(dict.lookupOrDefault<Foam::scalar>("Pv", 0)),       // Distal pressure
    rho_(dict.lookupOrDefault<Foam::scalar>("rho", 1)),     // Fluid density
    Qoo_(dict.lookupOrDefault<Foam::scalar>("Qoo", 0)),     // Historical flow rate values
    Qo_(dict.lookupOrDefault<Foam::scalar>("Qo", 0)),
    Qn_(dict.lookupOrDefault<Foam::scalar>("Qn", 0)),
    Poo_(dict.lookupOrDefault<Foam::scalar>("Poo", 0)),     // Historical pressure values
    Po_(dict.lookupOrDefault<Foam::scalar>("Po", 0)),
    Pn_(dict.lookupOrDefault<Foam::scalar>("Pn", 0)),
    Pioo_(dict.lookupOrDefault<Foam::scalar>("Pioo", 0)),   // Historical intermediate pressure values
    Pio_(dict.lookupOrDefault<Foam::scalar>("Pio", 0)),
    Pin_(dict.lookupOrDefault<Foam::scalar>("Pin", 0)),
    Pimoo_(dict.lookupOrDefault<Foam::scalar>("Pimoo", 0)), // Historical intramyocardial pressure
    Pimo_(dict.lookupOrDefault<Foam::scalar>("Pimo", 0)),
    Pimn_(dict.lookupOrDefault<Foam::scalar>("Pimn", 0)),
    timeIndex_(-1) // Initialize time index
{
    Info << "\n\nApplying coronaryOutletPressure BC on patch: " << patch().name() << endl;

    //- Extract differencing scheme from dictionary and validate
    word schemeStr = dict.lookupOrDefault<word>("diffScheme", "secondOrder");

    Info << "Differencing Scheme for Coronary LPN Model: " << schemeStr << endl;

    if (schemeStr == "firstOrder")
    {
        diffScheme_ = firstOrder; // First-order scheme
    }
    else if (schemeStr == "secondOrder")
    {
        diffScheme_ = secondOrder; // Second-order scheme
    }
    else
    {
        FatalErrorInFunction << "\n\nUnknown Differencing Scheme: " << schemeStr
                             << "\nValid Differencing Schemes (diffScheme) are : \n\n"
                             << " 2 \n ( \n firstOrder \n secondOrder \n ) \n"
                             << exit(FatalError);
    }

    //- Extract file name for intramyocardial pressure (Pim) data
    word PimFileStr = dict.lookupOrDefault<fileName>("PimFile", "PimData");

    PimFile_ = PimFileStr;

    //- Load the Pim data from the specified file
    loadPimData(PimFile_);

    //- Print coronary LPN model properties for debugging/verification
    Info << "\n\nProperties of Coronary LPN Model: \n"
            << "Differencig Scheme: " << diffScheme_ << " \n"
            << "Arterial Resistance: Ra = " << Ra_ << " kgm^-4s^-1 \n"
            << "Micro-arterial Resistance: Ram = " << Ram_ << " kgm^-4s^-1 \n"
            << "Veinous Resistance: Rv = " << Rv_ << " kgm^-4s^-1 \n"
            << "Micro-veinous Resistance: Rvm = " << Rvm_ << " kgm^-4s^-1 \n"
            << "Arterial Compliance: Ca = " << Ca_ << " m^4s^2kg^-1 \n"
            << "Intramyocardial Compliance: Cim = " << Cim_ << " m^4s^2kg^-1 \n"
            << "Path to Pim Data File: " << PimFile_ << " \n"
            << "Pim Scaling Factor: PimScaling: " << PimScaling_ << " \n"
            << "Pressure Distal to Rv: Pv = " << Pv_ << " Pa \n"
            << "Density: rho = " << rho_ << " kgm^-3 \n" << endl;
}

//- Constructor: Map an existing field onto a new patch
Foam::coronaryOutletPressureFvPatchScalarField::
coronaryOutletPressureFvPatchScalarField
(
    const coronaryOutletPressureFvPatchScalarField& ptf, // Existing field to map from
    const fvPatch& p,                                   // New patch
    const DimensionedField<Foam::scalar, volMesh>& iF,  // Internal field
    const fvPatchFieldMapper& mapper                    // Mapper for field data
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),    // Call base class mapping constructor
    Ra_(ptf.Ra_),                                       // Copy resistance parameter Ra
    Ram_(ptf.Ram_),                                     // Copy modified resistance Ram
    Rv_(ptf.Rv_),                                       // Copy venous resistance Rv
    Rvm_(ptf.Rvm_),                                     // Copy modified venous resistance Rvm
    Ca_(ptf.Ca_),                                       // Copy arterial compliance Ca
    Cim_(ptf.Cim_),                                     // Copy impedance compliance Cim
    PimScaling_(ptf.PimScaling_),                       // Copy scaling for Pim values
    Pv_(ptf.Pv_),                                       // Copy distal pressure Pv
    rho_(ptf.rho_),                                     // Copy density rho
    Qoo_(ptf.Qoo_), Qo_(ptf.Qo_), Qn_(ptf.Qn_),         // Copy flow rates Qoo, Qo, Qn
    Poo_(ptf.Poo_), Po_(ptf.Po_), Pn_(ptf.Pn_),         // Copy pressures Poo, Po, Pn
    Pioo_(ptf.Pioo_), Pio_(ptf.Pio_), Pin_(ptf.Pin_),   // Copy intermediate pressures Pioo, Pio, Pin
    Pimoo_(ptf.Pimoo_), Pimo_(ptf.Pimo_), Pimn_(ptf.Pimn_), // Copy Pim intermediate pressures
    diffScheme_(ptf.diffScheme_),                       // Copy differencing scheme
    PimFile_(ptf.PimFile_),                             // Copy Pim data file path
    PimData_(ptf.PimData_),                             // Copy Pim data storage
    timeIndex_(ptf.timeIndex_)                         // Copy current time index
{}

//- Copy constructor: Create a copy of an existing field
Foam::coronaryOutletPressureFvPatchScalarField::
coronaryOutletPressureFvPatchScalarField
(
    const coronaryOutletPressureFvPatchScalarField& copsf // Field to copy
)
:
    fixedValueFvPatchScalarField(copsf),               // Call base class copy constructor
    Ra_(copsf.Ra_),                                    // Copy resistance parameter Ra
    Ram_(copsf.Ram_),                                  // Copy modified resistance Ram
    Rv_(copsf.Rv_),                                    // Copy venous resistance Rv
    Rvm_(copsf.Rvm_),                                  // Copy modified venous resistance Rvm
    Ca_(copsf.Ca_),                                    // Copy arterial compliance Ca
    Cim_(copsf.Cim_),                                  // Copy impedance compliance Cim
    PimScaling_(copsf.PimScaling_),                    // Copy scaling for Pim values
    Pv_(copsf.Pv_),                                    // Copy distal pressure Pv
    rho_(copsf.rho_),                                  // Copy density rho
    Qoo_(copsf.Qoo_), Qo_(copsf.Qo_), Qn_(copsf.Qn_),  // Copy flow rates Qoo, Qo, Qn
    Poo_(copsf.Poo_), Po_(copsf.Po_), Pn_(copsf.Pn_),  // Copy pressures Poo, Po, Pn
    Pioo_(copsf.Pioo_), Pio_(copsf.Pio_), Pin_(copsf.Pin_), // Copy intermediate pressures
    Pimoo_(copsf.Pimoo_), Pimo_(copsf.Pimo_), Pimn_(copsf.Pimn_), // Copy Pim intermediate pressures
    diffScheme_(copsf.diffScheme_),                     // Copy differencing scheme
    PimFile_(copsf.PimFile_),                           // Copy Pim data file path
    PimData_(copsf.PimData_),                           // Copy Pim data storage
    timeIndex_(copsf.timeIndex_)                        // Copy current time index
{}

//- Copy constructor with new internal field reference
Foam::coronaryOutletPressureFvPatchScalarField::
coronaryOutletPressureFvPatchScalarField
(
    const coronaryOutletPressureFvPatchScalarField& copsf, // Field to copy
    const DimensionedField<scalar, volMesh>& iF            // New internal field reference
)
:
    fixedValueFvPatchScalarField(copsf, iF),           // Call base class copy with new field
    Ra_(copsf.Ra_),                                    // Copy resistance parameter Ra
    Ram_(copsf.Ram_),                                  // Copy modified resistance Ram
    Rv_(copsf.Rv_),                                    // Copy venous resistance Rv
    Rvm_(copsf.Rvm_),                                  // Copy modified venous resistance Rvm
    Ca_(copsf.Ca_),                                    // Copy arterial compliance Ca
    Cim_(copsf.Cim_),                                  // Copy impedance compliance Cim
    PimScaling_(copsf.PimScaling_),                    // Copy scaling for Pim values
    Pv_(copsf.Pv_),                                    // Copy distal pressure Pv
    rho_(copsf.rho_),                                  // Copy density rho
    Qoo_(copsf.Qoo_), Qo_(copsf.Qo_), Qn_(copsf.Qn_),  // Copy flow rates Qoo, Qo, Qn
    Poo_(copsf.Poo_), Po_(copsf.Po_), Pn_(copsf.Pn_),  // Copy pressures Poo, Po, Pn
    Pioo_(copsf.Pioo_), Pio_(copsf.Pio_), Pin_(copsf.Pin_), // Copy intermediate pressures
    Pimoo_(copsf.Pimoo_), Pimo_(copsf.Pimo_), Pimn_(copsf.Pimn_), // Copy Pim intermediate pressures
    diffScheme_(copsf.diffScheme_),                    // Copy differencing scheme
    PimFile_(copsf.PimFile_),                          // Copy Pim data file path
    PimData_(copsf.PimData_),                          // Copy Pim data storage
    timeIndex_(copsf.timeIndex_)                      // Copy current time index
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Update coefficients for boundary condition
void Foam::coronaryOutletPressureFvPatchScalarField::updateCoeffs()
{
    //- Check if coefficients are already updated
    if (updated())
    {
        return;
    }

    //- Get time step size and current simulation time
    const scalar dt = db().time().deltaTValue(); // Time step size
    const scalar currentTime = db().time().value(); // Current time

    //- Retrieve the flux field (phi) from the database
    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>("phi");

    //- Retrieve the boundary pressure field from the database
    const fvPatchField<scalar>& p = 
        db().lookupObject<volScalarField>("p").boundaryField()[this->patch().index()];

    //- Calculate total patch area
    scalar area = gSum(patch().magSf());

    //- Update pressure and flow rate histories if time index has advanced
    if (db().time().timeIndex() != timeIndex_)
    {
        timeIndex_ = db().time().timeIndex(); // Update time index

        //- Update flow rate history
        Qoo_ = Qo_;
        Qo_ = Qn_;
        Qn_ = gSum(phi.boundaryField()[patch().index()]); // Compute current flow rate

        //- Update pressure history
        Poo_ = Po_;
        Po_ = Pn_;
        Pn_ = gSum(p * patch().magSf()) / area; // Compute current pressure (area-weighted average)

        //- Update Pim history
        Pimoo_ = Pimo_;
        Pimo_ = Pimn_;
        Pimn_ = interpolatePim(currentTime); // Interpolate to compute current Pim

        //- Update intermediate pressure history
        Pioo_ = Pio_;
        Pio_ = Pin_;
    }

    //- Perform calculations only if time index exceeds 1
    if (db().time().timeIndex() > 1)
    {
        //- Choose differencing scheme
        switch (diffScheme_)
        {
            case firstOrder: // First-order scheme

                //- Compute intermediate pressure at current time
                Pin_ = (Rvm_ + Rv_) * Qn_ 
                    - (Rvm_ + Rv_) * Ca_ * ((Pn_ - Po_) / dt)
                    + (Rvm_ + Rv_) * Cim_ * ((Pimn_ - Pimo_) / dt)
                    + Pv_
                    + (Rvm_ + Rv_) * Cim_ * Pio_ / dt;

                Pin_ /= (1.0 + (Rvm_ + Rv_) * Cim_ / dt) + SMALL;

                //- Compute boundary pressure at current time
                Pn_ = (Ram_ + Ra_) * Qn_
                    + Ram_ * Ra_ * Ca_ * ((Qn_ - Qo_) / dt)
                    + Pin_ 
                    + Ram_ * Ca_ * (Po_ / dt);

                Pn_ /= (1.0 + Ram_ * Ca_ / dt) + SMALL;

                break;

            case secondOrder: // Second-order scheme
                
                //- Compute intermediate pressure at current time
                Pin_ = (Rvm_ + Rv_) * Qn_ 
                    - (Rvm_ + Rv_) * Ca_ * ((3 * Pn_ - 4 * Po_ + Poo_) / (2 * dt))
                    + (Rvm_ + Rv_) * Cim_ * ((3 * Pimn_ - 4 * Pimo_ + Pimoo_) / (2 * dt))
                    + Pv_
                    - (Rvm_ + Rv_) * Cim_ * (Pioo_ - 4 * Pio_) / (2 * dt);

                Pin_ /= (1.0 + 3 * (Rvm_ + Rv_) * Cim_ /(2 * dt)) + SMALL;

                //- Compute boundary pressure at current time
                Pn_ = (Ram_ + Ra_) * Qn_
                    + Ram_ * Ra_ * Ca_ * ((3 * Qn_ - 4 * Qo_ + Qoo_) / (2 * dt))
                    + Pin_ 
                    - Ram_ * Ca_ * (Poo_ - 4 * Po_) / (2 * dt);

                Pn_ /= (1.0 + 3 * Ram_ * Ca_ / (2 * dt)) + SMALL;

                break;

            default:
                //- Handle unknown schemes
                FatalErrorInFunction << "Unknown differencing scheme!" << exit(FatalError);
        }
    }

    //- Apply implicit update to the boundary field
    operator==(Pn_ / (rho_ + SMALL));

    //- Call base class function to finalize the update
    fixedValueFvPatchScalarField::updateCoeffs();
}

// Helper methods

//- Interpolates the intramyocardial pressure (Pim) at a given time
Foam::scalar Foam::coronaryOutletPressureFvPatchScalarField::interpolatePim(const scalar& time) const
{
    //- Ensure PimData_ is populated
    if (PimData_.empty())
    {
        FatalErrorInFunction << "PimData_ is empty. Cannot interpolate!" << exit(FatalError);
    }

    //- Get the start and end times of the Pim data
    const scalar tStart = PimData_.front().first; // Start time of the dataset
    const scalar tEnd = PimData_.back().first;   // End time of the dataset

    //- Calculate the effective time using modulo to handle periodicity
    scalar effectiveTime = tStart + std::fmod(time - tStart, tEnd - tStart);
    if (effectiveTime < tStart)
    {
        effectiveTime += (tEnd - tStart); // Adjust for negative modulo results
    }

    //- Perform linear interpolation to find the corresponding Pim value
    for (size_t i = 1; i < PimData_.size(); ++i)
    {
        if (PimData_[i - 1].first <= effectiveTime && PimData_[i].first >= effectiveTime)
        {
            scalar t1 = PimData_[i - 1].first;      // Previous time point
            scalar t2 = PimData_[i].first;          // Current time point
            scalar Pim1 = PimData_[i - 1].second;   // Pim value at previous time
            scalar Pim2 = PimData_[i].second;       // Pim value at current time

            //- Linear interpolation formula
            return PimScaling_ * (Pim1 + (Pim2 - Pim1) * (effectiveTime - t1) / (t2 - t1));
        }
    }

    //- Fallback case: return the last Pim value (should not typically be reached)
    return PimData_.back().second;
}

//- Loads intramyocardial pressure data (PimData_) from a file
void Foam::coronaryOutletPressureFvPatchScalarField::loadPimData(const word& fileName)
{
    //- Clear existing data
    PimData_.clear();

    //- Expand environment variables in the file name
    std::string expandedFileName = expandEnvironmentVariables(fileName);

    //- Open the file
    std::ifstream file(expandedFileName.c_str());

    //- Check if the file was successfully opened
    if (!file.is_open())
    {
        FatalErrorInFunction << "Cannot open intramyocardial pressure file: " << fileName << exit(FatalError);
    }

    std::string line;
    while (std::getline(file, line)) // Read the file line by line
    {
        //- Remove parentheses from the line
        line.erase(std::remove(line.begin(), line.end(), '('), line.end());
        line.erase(std::remove(line.begin(), line.end(), ')'), line.end());

        //- Skip empty or invalid lines
        if (line.empty()) continue;

        //- Parse the line into time and Pim values
        std::istringstream iss(line);
        scalar time, Pim;
        if (iss >> time >> Pim) // If parsing is successful
        {
            PimData_.emplace_back(time, Pim); // Store the time and Pim values
        }
    }

    file.close(); // Close the file
}

//- Expands environment variables in the given string
std::string Foam::coronaryOutletPressureFvPatchScalarField::expandEnvironmentVariables(const std::string& input)
{
    std::string result = input; // Start with the input string
    size_t pos = 0;

    //- Look for occurrences of "$" followed by an alphabetic character (environment variable pattern)
    while ((pos = result.find("$", pos)) != std::string::npos)
    {
        size_t endPos = result.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789", pos + 1);
        if (endPos == std::string::npos)
        {
            //- No valid environment variable found after "$"
            break;
        }

        //- Extract the environment variable name
        std::string varName = result.substr(pos + 1, endPos - pos - 1);

        //- Get the environment variable value from the system
        const char* envValue = std::getenv(varName.c_str());

        if (envValue) // If the environment variable exists
        {
            //- Replace the variable with its value
            result.replace(pos, endPos - pos, envValue);
            pos += std::string(envValue).size(); // Move past the replaced value
        }
        else
        {
            //- If the environment variable is not found, skip to the next "$"
            pos = endPos;
        }
    }
    return result; // Return the string with expanded variables
}

//- Write the boundary condition properties to an output stream
void Foam::coronaryOutletPressureFvPatchScalarField::write(Ostream& os) const
{
    //- Call base class method to write common properties
    fvPatchScalarField::write(os);

    //- Write each parameter with its keyword and value
    os.writeKeyword("Ra") << Ra_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ram") << Ram_ << token::END_STATEMENT << nl;
    os.writeKeyword("Rv") << Rv_ << token::END_STATEMENT << nl;
    os.writeKeyword("Rvm") << Rvm_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ca") << Ca_ << token::END_STATEMENT << nl;
    os.writeKeyword("Cim") << Cim_ << token::END_STATEMENT << nl;
    os.writeKeyword("PimFile") << PimFile_ << token::END_STATEMENT << nl;
    os.writeKeyword("PimScaling") << PimScaling_ << token::END_STATEMENT << nl;
    os.writeKeyword("Pv") << Pv_ << token::END_STATEMENT << nl;
    os.writeKeyword("rho") << rho_ << token::END_STATEMENT << nl;

    //- Write the differencing scheme
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
    os.writeKeyword("Qoo") << Qoo_ << token::END_STATEMENT << nl;
    os.writeKeyword("Qo") << Qo_ << token::END_STATEMENT << nl;
    os.writeKeyword("Qn") << Qn_ << token::END_STATEMENT << nl;

    //- Write pressure history parameters
    os.writeKeyword("Pimoo") << Pimoo_ << token::END_STATEMENT << nl;
    os.writeKeyword("Pimo") << Pimo_ << token::END_STATEMENT << nl;
    os.writeKeyword("Pimn") << Pimn_ << token::END_STATEMENT << nl;
    os.writeKeyword("Pioo") << Pioo_ << token::END_STATEMENT << nl;
    os.writeKeyword("Pio") << Pio_ << token::END_STATEMENT << nl;
    os.writeKeyword("Pin") << Pin_ << token::END_STATEMENT << nl;
    os.writeKeyword("Poo") << Poo_ << token::END_STATEMENT << nl;
    os.writeKeyword("Po") << Po_ << token::END_STATEMENT << nl;
    os.writeKeyword("Pn") << Pn_ << token::END_STATEMENT << nl;

    //- Write the boundary field value entry to the output stream
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Registration
namespace Foam
{
    //- Define the boundary field type for this class
    makePatchTypeField
    (
        fvPatchScalarField,                       // Base class
        coronaryOutletPressureFvPatchScalarField  // Derived class
    );
}

// ************************************************************************* //