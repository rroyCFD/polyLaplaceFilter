## Polynomial Laplace filter (polyLaplace) for Regularization and LES

### Author:
- Rajib Roy
- University of Wyoming
- rroy@uwyo.edu, roy.rajib@live.com

### Description

Polynomial Laplace filter for advection term regularization and LES. The polynomial coefficients are derived following the publication by Trias et al.

    Trias, F. X., & Verstappen, R. W. C. P. (2011). 
    On the construction of discrete filters for symmetry-preserving regularization models. Computers and Fluids, 40(1), 139–148. 
    https://doi.org/10.1016/j.compfluid.2010.08.015

### Dictionary template

#### Regularization (fvSolution sub-dictionary)

	regularization
	{
	    regOrder        C6; // C4 C2
	    filter          polyLaplace;

	    // epsilon (filter to grid spacing ratio) 3
	    d1              0.375; // 3/8
	    d2              0.0375; // 3/80

	    // epsilon (filter to grid spacing ratio) 2
	    // d1          0.16666667; // 1/6
	    // d2          0.00416667; //1/240

	}

#### LES filter (turbulenceProperties sub-dictionary)

	LESModelCoeffs
    {
        filter      polyLaplace;
        d1          0.16666667;  // 0.375;
        d2          0.00416667;  //  0.0375;
    }

### Make

*Location:* $WM_PROJECT_USER_DIR/src/TurbulenceModels/turbulenceModels/LES/LESfilters/polyLaplaceFilter

#### files
	
	LESfilters = LES/LESfilters
	$(LESfilters)/polyLaplaceFilter/polyLaplaceFilter.C

	LIB = $(FOAM_USER_LIBBIN)/libturbulenceModels_$(USER)

#### options

	EXE_INC = \
	    -I$(LIB_SRC)/TurbulenceModels/turbulenceModels/lnInclude \
	    -I$(LIB_SRC)/TurbulenceModels/incompressible/lnInclude \
	    -I$(LIB_SRC)/transportModels \
	    -I$(LIB_SRC)/finiteVolume/lnInclude \
	    -I$(LIB_SRC)/meshTools/lnInclude \

	LIB_LIBS = \
	    -lincompressibleTransportModels \
	    -lturbulenceModels \
	    -lfiniteVolume \
	    -lmeshTools