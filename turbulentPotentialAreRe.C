/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "turbulentPotentialAreRe.H"
#include "addToRunTimeSelectionTable.H"
#include "backwardsCompatibilityWallFunctions.H"
#include "components.H"
#include "fvCFD.H"
#include "volFields.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Hello

namespace Foam
{ 
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(turbulentPotentialAreRe, 0);
addToRunTimeSelectionTable(RASModel, turbulentPotentialAreRe, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


tmp<volScalarField> turbulentPotentialAreRe::Ts() const
{ 
    if(tslimiter_ == "true")
    {
        return max(k_/(epsilon_ + epsilonSmall_), 6.0*sqrt(nu()/(epsilon_ + epsilonSmall_)));
    }
    
    return ((k_+k0_)/(epsilon_ + epsilonSmall_));
}


tmp<volScalarField> turbulentPotentialAreRe::Ls() const
{
    
    volScalarField trueL = pow(k_+k0_, 1.5)/(epsilon_ + epsilonSmall_);
    
    if(lslimiter_ == "true")
    {
        return cL1_*max(pow(k_+k0_, 1.5)/(epsilon_ + epsilonSmall_),cL2_*pow(pow3(nu())/(epsilon_ + epsilonSmall_),0.25));
    }
    
    //Info << "Max trueL: " << gMax(trueL) << " Min trueL: " << gMin(trueL) << endl;
    //Info << "Dims: " << trueL.dimensions() << endl;
    
    return pow(k_+k0_, 1.5)/(epsilon_ + epsilonSmall_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentPotentialAreRe::turbulentPotentialAreRe
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    RASModel(typeName, U, phi, lamTransportModel),


    cEp1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEp1",
            coeffDict_,
            1.45
        )
    ),
    cEp2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEp2",
            coeffDict_,
            1.83
        )
    ),
    cEpType_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEpType",
            coeffDict_,
            1.0
        )
    ),
    cEp3_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEp3",
            coeffDict_,
            0.15
        )
    ),
    cEp4_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEp4",
            coeffDict_,
            0.05
        )
    ),
    cPe_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cPe",
            coeffDict_,
            0.1
        )
    ),
    cD1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cD1",
       	    coeffDict_,
            0.5
        )
    ),
    cD2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cD2",
       	    coeffDict_,
            0.33
        )
    ),
    cD3_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cD3",
       	    coeffDict_,
            3.0
        )
    ),
    cD4_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cD4",
       	    coeffDict_,
            2.4
        )
    ),
    cPm_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cPm",
            coeffDict_,
            1.0
        )
    ),
    cP1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP1",
            coeffDict_,
            2.0
        )
    ),
    cP2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP2",
            coeffDict_,
            0.6
        )
    ),
    cP3_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP3",
            coeffDict_,
            0.42
        )
    ),
    cP4_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP4",
            coeffDict_,
            0.85714
        )
    ),
    cP5_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cP5",
            coeffDict_,
            0.0
        )
    ),
    cGn_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cGn",
            coeffDict_,
            0.05
        )
    ),
    cGw_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cGw",
            coeffDict_,
            3.0
        )
    ),
    cLn_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cLn",
            coeffDict_,
            0.333
        )
    ),
    cLw_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cLw",
            coeffDict_,
            0.0
        )
    ),
    cN1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cN1",
            coeffDict_,
            1.0
        )
    ),
    cL1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cL1",
            coeffDict_,
            0.36
        )
    ),
    cL2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cL2",
            coeffDict_,
            85.0
        )
    ),
    cMu_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cMu",
            coeffDict_,
            0.21
        )
    ),
    betaK_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "betaK",
            coeffDict_,
            0.09
        )
    ),
    cT_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cT",
            coeffDict_,
            0.02
        )
    ),
    cA_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cA",
            coeffDict_,
            1.0
        )
    ),
    cEhmM_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cEhmM",
            coeffDict_,
            10.0
        )
    ),
	cNF_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cNF",
            coeffDict_,
            1.0
        )
    ),
	pMix_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "pMix",
            coeffDict_,
            1.0
        )
    ),
	cPrK_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cPrK",
            coeffDict_,
            0.6
        )
    ),
	cPrP_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cPrP",
            coeffDict_,
            1.0
        )
    ),
	cNL_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cNL",
            coeffDict_,
            0.00001
        )
    ),
	nutRatMax_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "nutRatMax",
            coeffDict_,
            1.0e5
        )
    ),
    sigmaK_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaK",
            coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            0.833
        )
    ),
    sigmaPhi_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaPhi",
            coeffDict_,
            0.33
        )
    ),
    sigmaPsi_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaPsi",
            coeffDict_,
            1.0
        )
    ),
    sCk_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sCk",
            coeffDict_,
            0.0
        )
    ),
    sCeps_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sCeps",
            coeffDict_,
            0.0
        )
    ),
    sCphi_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sCphi",
            coeffDict_,
            0.0
        )
    ),
    sCpsi_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sCpsi",
            coeffDict_,
            0.0
        )
    ),
    sCs_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sCs",
            coeffDict_,
            0.33
        )
    ),
	prodType_
	(
		dimensionedScalar::lookupOrAddToDict
		(
			"prodType",
			coeffDict_,
			1.0
		)
	),
	dampType_
	(
		dimensionedScalar::lookupOrAddToDict
		(
			"dampType",
			coeffDict_,
			1.0
		)
	),
	gType_
	(
		dimensionedScalar::lookupOrAddToDict
		(
			"gType",
			coeffDict_,
			1.0
		)
	),
    initSqrt_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "initSqrt",
            coeffDict_,
            0.0
        )
    ),

   eqncEp1_
   (
       coeffDict_.lookup("eqncEp1")
   ),
   
   eqncEp2_
   (
       coeffDict_.lookup("eqncEp2")
   ),

   eqnEpsHat_
   (
       coeffDict_.lookup("eqnEpsHat")
   ),

   debugWrite_
   (
       coeffDict_.lookup("debugWrite")
   ),
   tslimiter_
   (
       coeffDict_.lookup("tslimiter")
   ),
   lslimiter_
   (
       coeffDict_.lookup("lslimiter")
   ),
   eqncMu_
   (
       coeffDict_.lookup("eqncMu")
   ),
   nutReal_
   (
       coeffDict_.lookup("nutReal")
   ),
   y_
   (
   mesh_
   ),
    
	k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    
	gradk_
    (
        IOobject
        (
            "gradk",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (fvc::grad(k_))
    ),
    
	epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    
	nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    
	nutNorm_
    (
        IOobject
        (
            "nutNorm",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (nut_/max(nut_))
    ),
    
	tpphi_
    (
        IOobject
        (
            "tpphi",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    
	tpphiSqrt_
    (
        IOobject
        (
            "tpphiSqrt",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
	
	vorticity_
    (
        IOobject
        (
            "vorticity",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (fvc::curl(U_))
    ),
    
	tppsi_
    (
        IOobject
        (
            "tppsi",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    phiActual_
    (
        IOobject
        (
            "phiActual",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (sqr(tpphiSqrt_)*k_)
    ),

    psiActual_
    (
        IOobject
        (
            "psiActual",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (tppsi_*k_)
    ),
    
	uGrad_
    (
        IOobject
        (
            "uGrad",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (fvc::grad(U_))
    ),
	
	epsHat_
    (
        IOobject
        (
            "epsHat",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (epsilon_/(k_ + k0_))
    ),
    
	kSqrt_
    (
        IOobject
        (
            "kSqrt",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (sqrt(k_))
    ),
	
	alpha_
    (
        IOobject
        (
            "alpha",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (1.0/(1.0 + 1.5*sqr(tpphiSqrt_)))
    ),

    gamma_
    (
        IOobject
        (
            "gamma",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (alpha_)
    ),

    lambda_
    (
        IOobject
        (
            "lambda",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (alpha_)
    ),
	
	pod_
    (
        IOobject
        (
            "pod",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        tpphi_
    ),

    cMuRe_
    (
        IOobject
        (
            "cMuRe",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        tpphi_
    ),
	
	phiSqrt_
    (
        IOobject
        (
            "phiSqrt",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (sqrt(tpphi_*k_))
    ),
    
	gradkSqrt_
    (
        IOobject
        (
            "gradkSqrt",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (fvc::grad(kSqrt_))
    ),

    cEp1eqn_
    (
        IOobject
        (
            "cEp1eqn",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (cEp1_ + cEp4_*(2.0*alpha_ - 1.0))
    ),

	cEp2eqn_
    (
        IOobject
        (
            "cEp2eqn",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (cEp2_ - 0.16*exp(-0.1*sqr(k_)/(nu()*epsilon_)))
    ),
    
	tpProd_
    (
        IOobject
        (
            "tpProd",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        ((2*nut_*magSqr(symm(fvc::grad(U_)))/k_))
    ),
    
	cP1eqn_
    (
        IOobject
        (
            "cP1eqn",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (2.0*(0.5+0.5*((tpProd_*k_)/epsilon_)))
    ),
    
	gradTpphi_
    (
        IOobject
        (
            "gradTpphi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (fvc::grad(sqr(tpphiSqrt_)))
    ),
    
	gradTppsi_
    (
        IOobject
        (
            "gradTppsi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (fvc::grad(tppsi_))
    ),
    
	tpProdSqr_
    (
        IOobject
        (
            "tpProdSqr",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (sqr(tppsi_ & vorticity_))
    ),
    
	tpProd3d_
    (
        IOobject
        (
            "tpProd3d",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (mag(psiActual_ ^ vorticity_))
    ),
    sigFk_
    (
        IOobject
        (
            "sigFk",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tpphi_
    ),
    sigFeps_
    (
        IOobject
        (
            "sigFeps",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tpphi_
    ),
    sigFphi_
    (
        IOobject
        (
            "sigFphi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tpphi_
    ),
    sigFpsi_
    (
        IOobject
        (
            "sigFpsi",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tpphi_
    )
{


	
	
    //*************************************//	
    // Eddy viscosity - diffusion only
    //*************************************//
	nut_ = cMu_*k_*k_*tpphi_/epsilon_;	
    nut_ = min(nut_,nutRatMax_*nu());        
    nut_.correctBoundaryConditions();
    bound(nut_,dimensionedScalar("minNut", nut_.dimensions(), SMALL));       

	
    //*************************************//	
    // Epsilon-hat
    //*************************************//   
    epsHat_ = epsilon_/(k_ + (cEhmM_*nu()*mag(gradkSqrt_)));
    bound(epsHat_,dimensionedScalar("minEpsHat", epsHat_.dimensions(), SMALL));
	 

    // Take sqrt of initial tpphi field
    //if(initSqrt_.value()==1.0){
    //    Info<< "Taking sqrt of tpphi" << endl;
    //    tpphi_ = sqrt(tpphi_ + SMALL);
    //   tpphi_.correctBoundaryConditions();
    //}


    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Not used but necessary for RAS Model
tmp<volSymmTensorField> turbulentPotentialAreRe::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_))
        )
    );
}

// Not used but necessary for RAS Model
tmp<volSymmTensorField> turbulentPotentialAreRe::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}

// Term that is directly added to the momentum equation
tmp<fvVectorMatrix> turbulentPotentialAreRe::divDevReff() const
{
    return
    (
       fvc::grad(phiActual_)
     + fvc::curl(psiActual_)
     + fvc::laplacian(nut_, U_, "laplacian(nut,U)")
     - fvm::laplacian(nuEff(), U_)
    );
}


bool turbulentPotentialAreRe::read()
{
    if (RASModel::read())
    {
        cEp1_.readIfPresent(coeffDict());
        cEp2_.readIfPresent(coeffDict());
        cEpType_.readIfPresent(coeffDict());
        cEp3_.readIfPresent(coeffDict());
        cEp4_.readIfPresent(coeffDict());
        cPe_.readIfPresent(coeffDict());
        cP1_.readIfPresent(coeffDict());
        cP2_.readIfPresent(coeffDict());
        cP3_.readIfPresent(coeffDict());
		cP4_.readIfPresent(coeffDict());
        cP5_.readIfPresent(coeffDict());
        cGn_.readIfPresent(coeffDict());
        cGw_.readIfPresent(coeffDict());
        cLn_.readIfPresent(coeffDict());
        cLw_.readIfPresent(coeffDict());
        cN1_.readIfPresent(coeffDict());
		cL1_.readIfPresent(coeffDict());
		cL2_.readIfPresent(coeffDict());
        cMu_.readIfPresent(coeffDict());
		cEhmM_.readIfPresent(coeffDict());
		cPrK_.readIfPresent(coeffDict());
		cPrP_.readIfPresent(coeffDict());
		cD1_.readIfPresent(coeffDict());
		cD2_.readIfPresent(coeffDict());
		cD3_.readIfPresent(coeffDict());
		cD4_.readIfPresent(coeffDict());
		cT_.readIfPresent(coeffDict());
		cA_.readIfPresent(coeffDict());
		cNL_.readIfPresent(coeffDict());
		cNF_.readIfPresent(coeffDict());
		sigmaK_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        sigmaPhi_.readIfPresent(coeffDict());
		sigmaPsi_.readIfPresent(coeffDict());
        sCk_.readIfPresent(coeffDict());
        sCeps_.readIfPresent(coeffDict());
        sCphi_.readIfPresent(coeffDict());
        sCpsi_.readIfPresent(coeffDict());
		prodType_.readIfPresent(coeffDict());
		gType_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void turbulentPotentialAreRe::correct()
{

    //**********************************************//	
    // Bounding values not already defined by model
    //**********************************************//
	
    const dimensionedScalar eH0("minEpsHat", epsHat_.dimensions(), ROOTVSMALL);
	const dimensionedScalar nut0("minNut", nut_.dimensions(), ROOTVSMALL);
	const dimensionedScalar nutSmall("smallNut", nut_.dimensions(), SMALL);
	const dimensionedScalar tph0("minTpphi", tpphi_.dimensions(), ROOTVSMALL);
    const dimensionedScalar tphs0("minTpphi", tpphiSqrt_.dimensions(), ROOTVSMALL);
	const dimensionedScalar L0("lMin", dimensionSet(0,1,0,0,0,0,0), ROOTVSMALL);
    const dimensionedScalar v0("v0", vorticity_.dimensions(), SMALL);

    if (mesh_.changing())
    {
        y_.correct();
        bound(k_, k0_);
        bound(epsilon_, epsilonSmall_);
		bound(tpphi_,tph0);
        bound(tpphiSqrt_,tphs0);
		bound(nut_,nut0);
    }
	
	
    RASModel::correct();

	
    if (!turbulence_)  
    {
        return;
    } 
	
	
    //*************************************//	
    // Timestep - for use in elliptic switch
    //*************************************//
	
    dimensionedScalar cTime = U_.mesh().time().value();
	word sMMdebug = runTime_.controlDict().lookup("showMaxMin");


    //*************************************//	
    // Vorticity and Gradient
    //*************************************//
    
	vorticity_ = fvc::curl(U_);
	uGrad_ = fvc::grad(U_); 




	
	//*************************************//	
    // Misc Terms
    //*************************************//

	const volVectorField gradPhi_("gradPhi", fvc::grad(phiActual_));		
	const volScalarField gradgradPhi_("gradgradPhi", fvc::laplacian(DphiEff(),phiActual_));

    gradTpphi_ = fvc::grad(sqr(tpphiSqrt_));

	const volVectorField gradTpphiSqrt("gradTpphiSqrt",fvc::grad(tpphiSqrt_));

    phiSqrt_ = sqrt(tpphi_*k_ + k0_);
	const volVectorField gradPhiSqrt_("gradPhiSqrt",fvc::grad(phiSqrt_));	

    kSqrt_ = sqrt(mag(k_)+k0_);
    bound(kSqrt_,dimensionedScalar("minKsqrt", kSqrt_.dimensions(), sqrt(ROOTVSMALL)));

    gradk_ = fvc::grad(k_);
    gradkSqrt_ = fvc::grad(kSqrt_);

    gradTppsi_ = fvc::grad(tppsi_);


	
	
    //*************************************//	
    // K Production
    //*************************************//	

    const volScalarField S2(2*magSqr(symm(fvc::grad(U_))));
    volScalarField G("RASModel::G", nut_*S2);

	//tpProd_ = pMix_*(tppsi_ & vorticity_) + (1.0-pMix_)*G/(k_);
    tpProd_ = (tppsi_ & vorticity_);
	tpProdSqr_ = sqr(tpProd_);
	tpProd3d_ = mag(psiActual_ ^ vorticity_);	

	
    //*************************************//   
    // Epsilon-hat (Inverse time scale)
    //*************************************//
    
    epsHat_ = epsilon_/(k_ + (cEhmM_*nu()*mag(gradkSqrt_)) + k0_);
    bound(epsHat_,eH0);


	//*************************************//	
    // Update Alpha
    //*************************************//
     
	alpha_ = 1.0/(1.0 + 1.5*sqr(tpphiSqrt_));


    pod_ = (psiActual_ & psiActual_)/(cMu_*sqr(tpphiSqrt_)*sqr(k_) + sqr(k0_));
    bound(pod_,tph0);

    volScalarField gammaWall("gammaWall", cGw_*nu()*(gradTpphiSqrt & gradTpphiSqrt)*k_/epsilon_); 
    volScalarField lambdaWall("gammaWall", cLw_*nu()*(gradTpphiSqrt & gradTpphiSqrt)*k_/epsilon_); 
    volScalarField nutRe("nutRe", cMu_*sqr(tpphiSqrt_)*k_*k_/epsilon_);
    gamma_ = 1.0/(1.0 + cGn_*(nutRe/nu()) + gammaWall);
    lambda_ = 1.0/(1.0 + cLn_*sqrt(nutRe/nu()) + lambdaWall);

    volScalarField IIb("IIb", sqr(2.0*alpha_-1.0) + 2.0*(tppsi_ & tppsi_));  
    bound(IIb, SMALL);

    //*************************************//   
    // Diffusion Sigmas 
    //*************************************//

    sigFk_ = sCk_*(sCs_ + (sigmaK_-sCs_)*pod_) + (1.0-sCk_)*sigmaK_*(alpha_/alpha_);
    sigFeps_ = sCeps_*(sCs_ + (sigmaEps_-sCs_)*pod_) + (1.0-sCeps_)*sigmaEps_*(alpha_/alpha_);
    sigFphi_ = sCphi_*(sCs_ + (sigmaPhi_-sCs_)*pod_) + (1.0-sCphi_)*sigmaPhi_*(alpha_/alpha_);
    sigFpsi_ = sCpsi_*(sCs_ + (sigmaPsi_-sCs_)*pod_) + (1.0-sCpsi_)*sigmaPsi_*(alpha_/alpha_);

    //Info<< "Min sigFk: " << gMin(sigFk_) << " Min sigFk: " << gMax(sigFk_)  <<endl;
    //Info<< "Min sigFeps: " << gMin(sigFeps_) << " Min sigFeps: " << gMax(sigFeps_)  <<endl;
    //Info<< "Min sigFphi: " << gMin(sigFphi_) << " Min sigFphi: " << gMax(sigFphi_)  <<endl;
    //Info<< "Min sigFpsi: " << gMin(sigFpsi_) << " Min sigFpsi: " << gMax(sigFpsi_)  <<endl;




    //*************************************//
    //Dissipation equation
    //*************************************//

    cEp1eqn_ = cEp1_ + cEp4_*(2.0*alpha_-1.0);
    cEp2eqn_ = cEp2_ - 0.16*exp(-0.1*sqr(k_)/(nu()*(epsilon_ + epsilonSmall_)));	

    if(cEpType_.value() == 1.0){
        cEp1eqn_ = cEp1_ + cEp4_*(2.0*alpha_-1.0);
        cEp2eqn_ = cEp2_ + 0.0*(alpha_/alpha_);
    }

    if(cEpType_.value() == 2.0){
        cEp1eqn_ = cEp1_ + 0.0*(alpha_/alpha_);
        cEp2eqn_ = cEp2_ + 0.0*(alpha_/alpha_);
    }

    if(cEpType_.value() == 3.0){
        cEp1eqn_ = cEp1_ + 0.0*(alpha_/alpha_);
        cEp2eqn_ = cEp2_ - 0.16*exp(-0.1*sqr(k_)/(nu()*(epsilon_ + epsilonSmall_)));    
    }

    if(cEpType_.value() == 4.0){
        cEp1eqn_ = cEp1_ + 0.0*(alpha_/alpha_);
        cEp2eqn_ = cEp2_ - cEp4_*exp(-cPe_*phiActual_*k_/(nu()*(epsilon_ + epsilonSmall_)));
    }

    tmp<fvScalarMatrix> epsEqn   
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      + fvm::SuSp(-fvc::div(phi_), epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
       cEp1eqn_*tpProd_*k_*epsHat_
     - fvm::Sp(cEp2eqn_*epsHat_,epsilon_)
    );

    epsEqn().relax();
    solve(epsEqn);
    bound(epsilon_,epsilonSmall_);
	
	



	
    //*************************************//
    // Turbulent kinetic energy equation
    //*************************************//
    
    tmp<fvScalarMatrix> kEqn
    (

        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      + fvm::SuSp(-fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        tpProd_*k_
      - fvm::Sp(epsilon_/(k_+k0_),k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_,k0_);
  
	

	
	

	
    //*************************************//
    // Fix for separated flows 
    //*************************************// 	

	volScalarField cTexp("cTexp", cT_*sqrt((((nu()/100.0)+nut_)/nu())));
	
    volScalarField transPhi("transPhi", cTexp*cA_*((2.0/3.0) - sqr(tpphiSqrt_))*tpProd_);	
	volVectorField transPsi("transPsi", cTexp*((1.0 - alpha_)*vorticity_ - cA_*tppsi_*tpProd_));	
	
	
    volScalarField psDamp("psDamp", nut_/(nut_ + cNF_*nu()));
    volScalarField psDampOne("psDampOne", nut_/(nut_ + nu()));


    volScalarField cP1eqn("cP1eqn", 1.0 + cP1_*psDamp);
 

	//*************************************//
    // Phi/K^(1/2) equation 
    //*************************************//

    tmp<fvScalarMatrix> tpphiSqrtEqn
    ( 
        fvm::ddt(tpphiSqrt_)
      + fvm::div(phi_, tpphiSqrt_)
      + fvm::SuSp(-fvc::div(phi_), tpphiSqrt_)
      - fvm::laplacian(DphiEff(), tpphiSqrt_) 
      == 
	  // Pressure Strain Slow + Fast
	    0.5*cP1eqn*(2.0*alpha_-1.0)*epsHat_*tpphiSqrt_
	  + 0.5*cP2_*tpProd_*tpphiSqrt_

	  // Dissipation
	  - fvm::Sp(0.5*(2.0*alpha_-1.0)*epsHat_,tpphiSqrt_)
	  
      // From K
      - fvm::Sp(0.5*tpProd_,tpphiSqrt_)

	  // Transition
      + transPhi
    );  

    tpphiSqrtEqn().relax();
    solve(tpphiSqrtEqn);
    bound(tpphiSqrt_,tph0); 

	
    tpphi_ = sqr(tpphiSqrt_);
    tpphi_.correctBoundaryConditions();
    bound(tpphi_,tph0);   
 
  



    //*************************************//   
    // Psi Equation
    //*************************************//
    
    tmp<fvVectorMatrix> tppsiEqn
    (
        fvm::ddt(tppsi_)
      + fvm::div(phi_, tppsi_)
      + fvm::SuSp(-fvc::div(phi_), tppsi_)
      - fvm::laplacian(DpsiEff(), tppsi_)

      ==

	  // Production
	    tpphi_*vorticity_
		
	  // Slow Pressure Strain
      - cP1eqn*(1.0-alpha_)*epsHat_*tppsi_

	  // Fast Pressure Strain
	  - cP2_*tpphi_*vorticity_
	  + cP2_*tpProd_*tppsi_
	  - cP3_*2.0*alpha_*tpProd_*tppsi_
      + cP4_*(2.0*alpha_-1.0)*tpphi_*vorticity_

      // Dissipation
      + (1.0-alpha_)*epsHat_*tppsi_
	  //+ cP4_*(2.0*alpha_-1.0)*epsHat_*tppsi_

      // Recirc regions
      + cP5_*(k_*(1.0-alpha_)/(nu()/cLn_ + sqrt(nut_*nu() + sqr(nut0))))*tppsi_

	  // From K Equation
      - tpProd_*tppsi_

	  // Transition Term
      + transPsi
    );


    tppsiEqn().relax();
    solve(tppsiEqn);

	

	


	//*************************************//
    // Calculate eddy viscosity
    //*************************************//
   
	//nut_ = cMu_*k_*k_*sqr(tpphi_)/epsilon_;
	//nut_ = min(nut_,nutRatMax_*nu()); 
	//nut_.correctBoundaryConditions();
    //bound(nut_,nut0); 

    //volScalarField cMuEqn("cMuEqn", min((2.0*cMu_*alpha_/alpha_),(psiActual_ & psiActual_)*epsilon_/(sqr(tpphi_)*sqr(k_)*mag(tpProd_*k_) + sqr(k0_)*epsilonSmall_)));
    //volScalarField cMuEqn("cMuEqn", (psiActual_ & psiActual_)*epsilon_/(sqr(tpphi_)*sqr(k_)*mag(tpProd_*k_) + sqr(k0_)*epsilonSmall_));

    //volScalarField cMuRe("cMuRe", cN1_*0.7*cMuEqn + 0.15*(1-alpha_)*gamma_*cMuEqn*cN1_ + 0.06*(1.0 + 2.5*(1.0-cN1_)));
    //volScalarField cMuRe("cMuRe", 0.15*cN1_ + lambda_*cMuEqn*cN1_ + cMu_*(alpha_/alpha_)*(1.0-cN1_));
    //volScalarField cMuRe("cMuRe", (0.08 + 0.5*(1.1 + (1.0-alpha_)*gamma_)*cMuEqn)*cN1_ + cMu_*(alpha_/alpha_)*(1.0-cN1_));
    //cMuRe_ =  min((2.0*cMu_*alpha_/alpha_), 1.0/(1.0 + (1.0-cMu_)/(cMuEqn + tph0)));

    //cMuRe_ =  cN1_*(0.12 + 1.6*(1.0-(2.0/3.0)*IIb)*lambda_/(1.0 + (1.0-cMu_)/(cMuEqn + tph0))) + (1.0-cN1_)*cMu_*(alpha_/alpha_);
    //cMuRe_ = (0.08 + 0.5*(1.1 + (1.0-alpha_)*gamma_)*cMuEqn)*cN1_ + cMu_*(alpha_/alpha_)*(1.0-cN1_);
 
    //nut_ = cMuRe_*sqr(tpphi_)*k_*k_/epsilon_; 

    //nut_ = 1.0/((1.0-0.8*sqrt(alpha_))*epsilon_/(cMu_*phiActual_*k_ + k0_*k0_) + 0.8*sqrt(alpha_)*tpProd_*k_/((psiActual_ & psiActual_) + (k0_*k0_)) );
    volScalarField magProd("magProd", mag((tpProd_*k_)));
    bound(magProd, epsilonSmall_);


    //nut_ = 1.0/((gradkSqrt_ & gradkSqrt_)/magProd + magProd/((psiActual_ & psiActual_) + (k0_*k0_)));
    // nut_ = ((0.12*gamma_ + 0.17)*sqr(tpphi_)*cN1_ + cMu_*(sqr(tpphi_)*(1.0-cN1_)))*k_/epsHat_;

    cMuRe_ = cN1_*(0.5*(0.12 + 0.37*lambda_)*tpphi_ + 0.5*(tppsi_ & tppsi_)) + (1.0-cN1_)*cMu_*tpphi_;

    nut_ = cMuRe_*k_*k_/epsilon_;    
    nut_ = min(nut_,nutRatMax_*nu());  
    nut_.correctBoundaryConditions();
    bound(nut_,nut0);   
 


    //*************************************//
    // Variables to add to NS
    //*************************************//
    phiActual_ = sqr(tpphiSqrt_)*k_;
    psiActual_ = tppsi_*k_;


	
    //*************************************//   
    // Output some max values
    //*************************************//
	
	if(sMMdebug == "true")
	{    
	volScalarField meanUz("meanUz",U_.component(2));
	volScalarField uTauSquared((nu() + nut_)*vorticity_.component(2));
	volVectorField tpphiVort(sqr(tpphiSqrt_)*vorticity_);
    volScalarField G("G", tpProd_*k_);

	
	Info<< "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << endl; 
    Info<< "Max cMuRe: " << gMax(cMuRe_) << " Min cMuRe: " << gMin(cMuRe_)  <<endl;
	Info<< "Max Epsilon: " << gMax(epsilon_) << " Min Epsilon: " << gMin(epsilon_) << " Max G: " << gMax(G) <<endl;
	Info<< "Max nut: " << gMax(nut_) << " Max K: " << gMax(k_) << " Max Phi: " << gMax(phiActual_) <<endl;
    Info<< "Max Psi: " << gMax(psiActual_) << " Min Psi: " << gMin(psiActual_)  <<endl;
    Info<< "Max Phi/k: " << gMax(tpphi_) << " Min Phi/k: " << gMin(tpphi_)  <<endl;
    Info<< "Max Phi/k Sqrt: " << gMax(tpphiSqrt_) << " Min Phi/k: " << gMin(tpphiSqrt_)  <<endl;
    Info<< "Max uTauSquared: " << gMax(uTauSquared) <<endl; 
	Info<< "Max vorticity: " << gMax(vorticity_) << " Min vorticity: " << gMin(vorticity_) <<endl;
	Info<< "Max Uz: " << gMax(meanUz) << " Min Uz: " << gMin(meanUz) << endl;
	Info<< "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << endl; 
	}

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
