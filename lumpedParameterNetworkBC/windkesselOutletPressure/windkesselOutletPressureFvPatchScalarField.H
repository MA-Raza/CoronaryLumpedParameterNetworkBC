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
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::windkesselOutletPressureFvPatchScalarField

Group


Description
    A boundary condition that models outlet pressure using the Windkessel 
    framework, commonly used in hemodynamics simulations to replicate the 
    physiological behavior of arterial systems.

    This boundary condition can simulate different Windkessel models of 
    varying complexity, including:
    
    (1) Resistive (R): A simple resistive boundary where pressure drop is
        proportional to flow.
        
    (2) WK2 (2-element Windkessel): A model with a distal resistance (Rd)
        and capacitance (C), commonly used to represent the compliance of 
        arteries without considering proximal resistance.

    (3) WK3 (3-element Windkessel): The classic 3-element model, which 
        includes proximal resistance (Rp), distal resistance (Rd), and 
        capacitance (C), and is the most widely used model in hemodynamics.

    (4) WK4Series (4-element Windkessel with inductance in series): An 
        extension of the 3-element model that adds inductance or inertance (L)
        accounting for the inertial effects of the blood in series with the 
        proximal resistance (Rp), making it more closely match patient-specific
        scenarios.

    (5) WK4Parallel (4-element Windkessel with inductance in parallel): Another 
        extension of the 3-element model, where inductance or inertance (L) is 
        placed in parallel to the proximal resistance (Rp), capturing additional 
        dynamics related to flow fluctuations and inertial effects.
    
    A distal pressure term (Pd) is added to all the models after the distal 
    resistance (Rd), making them applicable to cases where distal resistance 
    (Rd) is not grounded to zero pressure.

    The Resistive (R) model is based on the following equation for the 
    relationship between pressure (P) and flow rate (Q):

        P(t) = Rd * Q(t) + Pd

    Where:
        - Rd: Distal resistance
        - Pd: Distal pressure
        - Q(t): Flow rate
        - P(t): Pressure

    The 2-element Windkessel (WK2) model is based on the following ordinary
    differential equation (ODE) for the relationship between pressure (P) and
    flow rate (Q):
    
        P(t) = Rd * Q(t)
          - Rd * C * dP/dt
          + Pd

    Where:
        - Rd: Distal resistance
        - C: Compliance
        - Pd: Distal pressure
        - Q(t): Flow rate
        - P(t): Pressure

    The 3-element Windkessel (WK3) model is based on the following ordinary
    differential equation (ODE) for the relationship between pressure (P) and
    flow rate (Q):

        P(t) = (Rp + Rd) * Q(t)
          + (C * Rd * Rp) * dQ/dt
          + Pd
          - (C * Rd) * dP/dt

    Where:
        - Rd: Distal resistance
        - Rp: Proximal resistance
        - C: Compliance
        - Pd: Distal pressure
        - Q(t): Flow rate
        - P(t): Pressure

    The 4-element Windkessel model with inductance in series (WK4Series) is
    based on the following ordinary differential equation (ODE) for the 
    relationship between pressure (P) and flow rate (Q):

        P(t) = (Rp + Rd) * Q(t)
          + (L + Rp * Rd * C) * dQ/dt
          + Rd * C * L * d^2Q/dt^2
          + Pd
          - Rd * C * dP/dt

    Where:
        - Rd: Distal resistance
        - Rp: Proximal resistance
        - C: Compliance
        - L: Inductance or inertance
        - Pd: Distal pressure
        - Q(t): Flow rate
        - P(t): Pressure

    The 4-element Windkessel model with inductance in parallel (WK4Parallel) is 
    based on the following ordinary differential equation (ODE) for the 
    relationship between pressure (P) and flow rate (Q):

        P(t) = Rd * Q(t)
            + L * (1 + Rd / Rp) * dQ/dt
            + Rd * C * L * d^2Q/dt^2
            + Pd
            - ((L + Rp * Rd * C) / Rp) * dP/dt
            - (Rd * C * L / Rp) * d^2P/dt^2

    Where:
        - Rd: Distal resistance
        - Rp: Proximal resistance
        - C: Compliance
        - L: Inductance or inertance
        - Pd: Distal pressure
        - Q(t): Flow rate
        - P(t): Pressure

    This boundary condition calculates the outlet pressure based on volumetric 
    flow rates and the Windkessel parameters using either second-order backward
    differencing scheme or a simpler first-order backward differencing approach.

    This ODE is discretized using either a second-order backward Euler 
    scheme or a first-order scheme. The time derivative terms (dP/dt, dQ/dt,
    d^2P/dt^2, and d^2Q/dt^2) are handled differently depending on the selected
    differencing scheme.

    For first-order backward differencing scheme:
    
        dx/dt = (xn - xo) / dt

    and

        d^2x/dt^2 = (xn - 2 * xo + xoo) / (dt^2)

    Whereas, for second-order backward differencing scheme:

        dx/dt = (3*xn - 4*xo + xoo) / (2 * dt)

    and

        d^2x/dt^2 = (2 * xn - 5 * xo + 4 * xoo - xooo) / (dt^2)

    The boundary condition computes the pressure (Pn) at each time 
    step based on historical values of pressure and flow rates, 
    ensuring the boundary accurately models physiological behavior.

Key Features
    - Supports both first-order and second-order backward differencing scheme for
      pressure updates.
    - Updates historical pressure and flow rate values to reflect Windkessel dynamics.
    - Includes parameters such as Rp, Rd, Rc, L, and Pd, which are customizable 
      through a dictionary.

Usage
    Depending on the model, specify this boundary condition in the boundary file (0/p) as
    follows:

        <patchName>
        {
        type                windkesselOutletPressure;
        windkesselModel     <Resistive, WK2, WK3, WK4Series or WK4Parallel>; // Type of model, Default: WK3
        Rp                  <value>;    // Proximal resistance in [kgm^-4s^-1], Default: 1
        Rd                  <value>;    // Distal resistance in [kgm^-4s^-1], Default: 1
        C                   <value>;    // Capacitance or compliance in [m^4s^2kg^-1], Default: 1
        L                   <value>;    // Inductance or inertance in [kgm^-4], Default: 1
        Pd                  <value>;    // Distal pressure in [Pa], Default: 0
        rho					<value>;    // Fluid density in [kgm^-3], Default: 1
        diffScheme          <firstOrder or secondOrder>; // Differencing scheme, Default: secondOrder
        value               uniform 0;
        }
    
Note
    1) For the Resistive (R) boundary condition set windkesselModel = Resistive, 
       and give values of Rd and Pd.
    2) For the 2-element Windkessel (WK2) model set windkesselModel = WK2, and give 
       values of Rd, C, and Pd. Choose diffScheme.
    3) For the 3-element Windkessel (WK3) model set windkesselModel = WK3, and give 
       values of Rp, Rd, C, and Pd. Choose diffScheme.
    4) For the 4-element Windkessel (WK4Series) model with series inductance set 
       windkesselModel = WK4Series, and give values of Rp, Rd, C, L, and Pd. Choose 
       diffScheme.
    5) For the 4-element Windkessel (WK4Parallel) model with parallel inductance set 
       windkesselModel = WK4Parallel, and give values of Rp, Rd, C, L, and Pd. Choose 
       diffScheme.

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
    windkesselOutletPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef windkesselOutletPressureFvPatchScalarField_H
#define windkesselOutletPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Enumerator for supported windkessel model types
enum WindkesselModel
{
    Resistive = 0,     // Simple resistive model (R)
    WK2 = 1,           // 2-element Windkessel model (RC)
    WK3 = 2,           // 3-element Windkessel model (RCR)
    WK4Series = 3,     // 4-element Windkessel with inductance in series (RCRL)
    WK4Parallel = 4    // 4-element Windkessel with inductance in parallel
};

//- Enumerator for for time differencing schemes
enum DifferencingScheme
{
    firstOrder = 0,    // First-order backward differencing scheme
    secondOrder = 1    // Second-order backward differencing scheme
};

/*---------------------------------------------------------------------------*\
            Class windkesselOutletPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

/**
 * @brief Boundary condition implementing hemodynamic outlet pressure 
 *        using windkessel models.
 *
 * This class models outlet pressure as a function of volumetric flow rate
 * using the windkessel framework. The framework includes models of varying
 * complexity:
 *
 * - Resistive (R)
 * - 2-element Windkessel (WK2 or RC)
 * - 3-element Windkessel (WK3 or RCR)
 * - 4-element Windkessel with inductance in series (WK4Series or RCRL-Series)
 * - 4-element Windkessel with inductance in parallel (WK4Parallel or RCRL-Parallel)
 *
 * The models incorporate parameters such as proximal resistance (Rp), distal
 * resistance (Rd), compliance (C), inductance (L), and distal pressure (Pd). 
 * The boundary condition supports explicit and implicit formulations for time 
 * stepping and allows customization of differencing schemes and parameter settings.
 */

class windkesselOutletPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

private:

    // Private data

        //- Windkessel parameters
        scalar Rp_;     // Proximal resistance (Rp)
        scalar Rd_;     // Distal resistance (Rd)
        scalar C_;      // Compliance or capacitance (C)
        scalar L_;      // Inductance (L)
        scalar Pd_;     // Distal pressure (Pd)
        scalar rho_;    // Fluid density (rho)

        //- Windkessel state variables
        scalar Pooo_;   // Pressure three time steps ago
        scalar Poo_;    // Pressure two time steps ago
        scalar Po_;     // Previous time-step pressure
        scalar Pn_;     // Current pressure
        scalar Qooo_;   // Flow rate three time steps ago
        scalar Qoo_;    // Flow rate two time steps ago
        scalar Qo_;     // Previous time-step flow rate
        scalar Qn_;     // Current flow rate

        //- Windkessel LPN configuration
        WindkesselModel windkesselModel_;   // Selected Windkessel model type
        DifferencingScheme diffScheme_;     // Selected time differencing scheme

        //- Time index
        label timeIndex_;   // Tracks the simulation time index

public:

    //- Runtime type information
    TypeName("windkesselOutletPressure");

    // Constructors

        //- Construct from patch and internal field
        windkesselOutletPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field, and dictionary
        windkesselOutletPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Map existing object onto a new patch
        windkesselOutletPressureFvPatchScalarField
        (
            const windkesselOutletPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        windkesselOutletPressureFvPatchScalarField
        (
            const windkesselOutletPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new windkesselOutletPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        windkesselOutletPressureFvPatchScalarField
        (
            const windkesselOutletPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new windkesselOutletPressureFvPatchScalarField(*this, iF)
            );
        }

    // Member functions

        //- Evaluation function to update coefficients for the boundary 
        //  condition based on flow rate and model parameters
        virtual void updateCoeffs();

        //- Write boundary condition data to output stream
        virtual void write(Ostream&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //