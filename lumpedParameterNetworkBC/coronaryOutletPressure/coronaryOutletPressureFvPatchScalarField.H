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

Class
    Foam::coronaryOutletPressureFvPatchScalarField

Group


Description
    A boundary condition that models outlet pressure using the lumped parameter
    network (LPN) model of downstream coronary vascular beds, specifically
    designed for simulating hemodynamics in coronary artery systems. The LPN 
    model is commonly used in cardiovascular simulations to capture the complex
    interaction between the arterial system and various lumped components representing
    resistance, compliance, and pressure dynamics.

    The model uses a system of ordinary differential equations (ODEs) to relate 
    pressure (P) and flow rate (Q) at the boundary, incorporating parameters 
    like arterial resistance (Ra), micro-arterial resistance (Ram), arterial compliance
    (Ca), intramyocardial compliance (Cim), intramyocardial pressure (Pim), veinous
    resistance (Rv) and micro-veinous resistance (Rvm). Moreover, a distal pressure term
    (Pv) is added distal to Rv on veinous side. The veinous compliance (Cv) is not considered
    to simplify the numerics. The LPN model used have been shown in the literature to produce
    the characteristic coronary flow and pressure curve efficiently.

    The system of ordinary differential equations (ODEs) relating pressure (P) and flow rate
    (Q) for the LPN model of downstream coronary vascular beds is given as follows:

        (1)
        Pi(t) = (Rvm + Rv) * Q(t) 
            - (Rvm + Rv) * Ca * dP/dt
            + (Rvm + Rv) * Cim * dPim/dt
            + Pv
            - (Rvm + Rv) * Cim * dPi/dt

        (2)
        P(t) = (Ram + Ra) * Q(t)
            + Ram * Ra * Ca * dQ/dt
            + Pi(t) 
            - Ram * Ca * dP/dt

    Where:
        - Ra:    Arterial resistance
        - Ram:   Micro-arterial resistance
        - Rv:    Veinous resistance
        - Rvm:   Micro-veinous resistance
        - Ca:    Arterial compliance
        - Cim:   Intramyocardial compliance
        - Pim:   Intramyocardial pressure
        - Pv:    Pressure distal to Rv
        - Q(t):  Flow rate
        - Pi(t): Intermediate pressure
        - P(t):  Coronary outlet pressure
    
    This boundary condition calculates the coronary outlet pressure based on 
    volumetric flow rates and LPN model parameters using either second-order
    backward differencing scheme or a simpler first-order backward differencing
    approach.

    These ODEs are discretized using either a second-order backward Euler 
    scheme or a first-order scheme. The time derivative terms (dP/dt, dQ/dt,
    dPi/dt, and dPim/dt) are handled differently depending on the selected
    differencing scheme.
    
    For first-order backward differencing scheme:
    
        dx/dt = (xn - xo) / dt

    Whereas, for second-order backward differencing scheme:

        dx/dt = (3*xn - 4*xo + xoo) / (2 * dt)

    The boundary condition computes the pressure (P) at each time 
    step based on historical values of pressure and flow rates, 
    ensuring the boundary accurately models physiological behavior.

Key Features
    - Supports both first-order and second-order backward differencing scheme for
      pressure updates.
    - Updates historical pressure and flow rate values to reflect LPN dynamics.
    - Includes parameters such as Ra, Ram, Ca, Cim, Rvm, Rv, and Pv which are 
      customizable through a dictionary.
    - The intramyocardial pressure time-series data (Pim(t)) can be provided through
      a file.
    - A dictionary entery (PimScaling) can be used to scale the provide Pim(t) data.
      This is particularly useful for the cases where either left or right ventricle
      pressure is used as Pim(t), and the multiplication with different scaling factors
      is required for left and right coronary arteries (LCA and RCA).

Usage

    Specify this boundary condition in the boundary file (0/p) as follows:

    <patchName>
    {
    type                coronaryOutletPressure;
    Ra                   <value>;   //Arterial resistance in [kgm^-4s^-1], Default: 1
    Ram                  <value>;   //Micro-arterial resistance in [kgm^-4s^-1], Default: 1
    Rv                   <value>;   //Veinous resistance in [kgm^-4s^-1], Default: 1
    Rvm                  <value>;   //Micro-venous resistance in [kgm^-4s^-1], Default: 0
    Ca                   <value>;   //Arterial compliance in [m^4s^2kg^-1], Default: 1
    Cim                  <value>;   //Intramyocardial compliance in in [m^4s^2kg^-1], Default: 1
    PimFile              "<Path to Pim Data File>"; //Path to Pim time-series data file, Default: PimData
    PimScaling           <value>;   // Intramyocardial pressure scaling factor, Default: 1
    Pv                   <value>;   //Pressure distal to Rv in [Pa], Default: 0
    rho					 <value>; 		 // Fluid density in [kgm^-3], Default: 1
    diffScheme           <firstOrder or secondOrder>; // Differencing scheme, Default: secondOrder
    value                uniform 0;
    }

Author
    Muhammad Ahmad Raza, University College Dublin, Ireland.

Cite as
    Raza, M. A.: Implementation of Lumped Parameter Network Boundary Conditions for the
    Patient-Specific CFD Simulations of Coronary Arteries in OpenFOAM. In Proceedings of CFD with
    OpenSource Software, 2024, Edited by Nilsson. H., http://dx.doi.org/10.17196/OS_CFD#YEAR_2024

See also
    Foam::fixedValueFvPatchField
    Foam::fvPatchField

SourceFiles
    coronaryOutletPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef coronaryOutletPressureFvPatchScalarField_H
#define coronaryOutletPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "scalar.H"
#include "word.H"
#include <vector>  // Required for std::vector
#include <utility> // Required for std::pair

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Enumerator for numerical differencing schemes used in the LPN model
enum DifferencingScheme
{
    firstOrder = 0, // First-order backward differencing scheme
    secondOrder = 1 // Second-order backward differencing scheme
};

/*---------------------------------------------------------------------------*\
            Class coronaryOutletPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

/**
 * @brief Custom boundary condition for simulating hemodynamic outlet flow 
 *        using a coronary lumped parameter network (LPN) model.
 *
 * This boundary condition models the outlet pressure based on a system of
 * ordinary differential equations (ODEs) capturing downstream coronary
 * vascular dynamics. It incorporates resistance, compliance, and 
 * time-varying intramyocardial pressure to simulate physiological behavior
 * of coronary artery systems.
 *
 * The implementation supports first-order and second-order backward
 * differencing schemes, customizable parameters, and external input for
 * intramyocardial pressure data.
 */

class coronaryOutletPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

private:

    // Private data

        //- Model parameters
        scalar Ra_;          // Arterial resistance
        scalar Ram_;         // Micro-arterial resistance
        scalar Rv_;          // Veinous resistance
        scalar Rvm_;         // Micro-veinous resistance
        scalar Ca_;          // Arterial compliance
        scalar Cim_;         // Intramyocardial compliance
        scalar PimScaling_;  // Scaling factor for intramyocardial pressure
        scalar Pv_;          // Pressure distal to Rv
        scalar rho_;         // Fluid density

        //- State variables for flow rates and pressures
        scalar Qoo_, Qo_, Qn_;       // Historical flow rate values
        scalar Poo_, Po_, Pn_;       // Historical pressure values
        scalar Pioo_, Pio_, Pin_;    // Historical intermediate pressure values
        scalar Pimoo_, Pimo_, Pimn_; // Historical intramyocardial pressure values

        //- LPN model configuration
        DifferencingScheme diffScheme_; // Numerical differencing scheme
        Foam::fileName PimFile_;      // File containing intramyocardial pressure data

        //- Intramyocardial pressure data
        std::vector<std::pair<scalar, scalar>> PimData_; // Time-series data for (time, Pim)

        //- Time index
        label timeIndex_;             // Current time index for updates

    // Helper methods

        //- Function to interpolate intramyocardial pressure (Pim) data at a given time
        scalar interpolatePim(const scalar& time) const;

        //- Function to load intramyocardial pressure (Pim) data from a specified file
        void loadPimData(const word& fileName);

        //- Function to expands environment variables in a given string
        static std::string expandEnvironmentVariables(const std::string& input);

public:

    //- Runtime type information
    TypeName("coronaryOutletPressure");

    // Constructors

        //- Construct from patch and internal field
        coronaryOutletPressureFvPatchScalarField
        (
            const fvPatch& patch,
            const DimensionedField<scalar, volMesh>& field
        );

        //- Construct from patch, internal field, and dictionary
        coronaryOutletPressureFvPatchScalarField
        (
            const fvPatch& patch,
            const DimensionedField<scalar, volMesh>& field,
            const dictionary& dict
        );

        //- Map existing object onto a new patch
        coronaryOutletPressureFvPatchScalarField
        (
            const coronaryOutletPressureFvPatchScalarField& other,
            const fvPatch& patch,
            const DimensionedField<scalar, volMesh>& field,
            const fvPatchFieldMapper& mapper
        );

        //- Copy constructor
        coronaryOutletPressureFvPatchScalarField
        (
            const coronaryOutletPressureFvPatchScalarField& other
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new coronaryOutletPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        coronaryOutletPressureFvPatchScalarField
        (
            const coronaryOutletPressureFvPatchScalarField& other,
            const DimensionedField<scalar, volMesh>& field
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new coronaryOutletPressureFvPatchScalarField(*this, iF)
            );
        }

    // Member functions

        //- Evaluation function to update coefficients for the boundary 
        //  condition based on flow rate and model parameters
        virtual void updateCoeffs();

        //- Write boundary condition data to output stream
        virtual void write(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

