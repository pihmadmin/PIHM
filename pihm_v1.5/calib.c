/*******************************************************************************
 * File        : calib.c                                                       *
 * Function    : defines calibration factors for physical parameters           *
 * Programmers : Yizhong Qu   @ Pennsylvania State Univeristy                  *
 *               Mukesh Kumar @ Pennsylvania State Univeristy                  *
 *               Gopal Bhatt  @ Pennsylvania State Univeristy                  *
 * Version     : 2.0 (July 10, 2007)                                           *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *                                                                             *
 * This code is free for users with research purpose only, if appropriate      *
 * citation is refered. However, there is no warranty in any format for this   *
 * product.                                                                    *
 *                                                                             *
 * For questions or comments, please contact the authors of the reference.     *
 * One who want to use it for other consideration may also contact Dr.Duffy    *
 * at cxd11@psu.edu.                                                           *
 *******************************************************************************/

//! @file calib.c Calibration parameters and function definitions

/* SUNDIALS Header Files */
#include "sundials_types.h"

#define FileName        "rhode";

#define satD_CALIB      0.0
#define br_CALIB        5.0
#define poros_CALIB     0.35        /* Multiplication factor for element porosity                              */
#define icsat_CALIB     0.45        /* Set the initial saturation of all the elements to this value            */
#define rivEle_CALIB    0.45        /* Set the initial saturation of elements near river segment to this value */

#define is_CALIB        20.0        /* Multiplicatoin factor for Maximum Interception Storage                  */
#define et0_CALIB       9.0         /* Multiplication factor for Rate of Evaporation from Canopy               */
#define mf_CALIB        1.50        /* Multiplication factor for Rate of Snow Melt                             */
#define tf_CALIB        1.0         /* Multiplication factor for Rate of Throughfall                           */

#define Vic_CALIB       0.5         /* Fractional factor for Rate of Infiltration Capacity                     */
#define rivK_CALIB      50.0        /* Multiplicatoin factor for Conductivity at River-Element Interface       */
#define Kh_CALIB        1.0         /* Multiplicatoin factor for Horizontal Conductivity of Elements           */
#define Rec_CALIB       1.0         /* Multiplication factor for Rate of groundwater recharge                  */
#define et2_CALIB       1.0         /* Multiplication factor for ET2                                           */
#define et1_CALIB       1.0         /* Multiplicatoin factor for ET1                                           */
#define sat_THRESH      0.685       /* Saturation Threshold above which macropores respond                     */
#define mp_MULTFH       3000000.0   /* Multiplicatoin factor for Horizontal Conductivity of Macropores         */
#define mp_MULTFV       3000000.0   /* Multiplication factor for Vertical Conductivity of Macropores           */
#define mpArea_CALIB    0.00        /* Fractional Area of macropore of any element                             */
#define ovl_THRESH_H    -1.0
#define ovl_THRESH_V    0.0
#define rzd_CALIB       0.20        /* Subtract value from Root Zone Depth                                     */

#define roughEle_CALIB  1.0         /* Multiplication factor for Manning's roughness of Element                */
#define roughRiv_CALIB  1.0         /* Multiplication factor for Manning's roughness of River Segments         */
#define rivCoeff_CALIB  1.0         /* Multiplication factor for River Segment Interpolation Coefficient       */
                                    /* (width incase of a rectangular river segment)                           */
#define rivDepth_CALIB  2.5         /* Multiplicatoin factor for Depth of a River segment                      */
#define alpha_CALIB     1.0         /* Multiplicatoin factor for Soil Alpha parameter                          */
#define set_MP          0.0         /* An Identifier to set if soil exibit macropore (1: Yes, 0: No)           */
#define lai_CALIB       1.00        /* Multiplication factor for Leaf Area Index                               */
#define vegfrac_CALIB   1.00        /* Multiplication factor for Vegitation fraction                           */
#define albedo_CALIB    1.25        /* Multiplication factor for Albedo                                        */

void setFileName(char *fileName)
//! Input File Name ID ("test" => test.mesh test.att ...)
{
	char tempFileName[100] = FileName;
	strcpy(fileName, tempFileName);
	return;
}

realtype setsatD_CALIB(){
	return satD_CALIB;
}
realtype setbr_CALIB(){
	return br_CALIB;
}
realtype setporos_CALIB()
//! Calibration factor for soil porosity parameter
{
	return poros_CALIB;
}
realtype seticsat_CALIB()
//! Inital saturation of elements
{
	return icsat_CALIB;
}
realtype setrivEle_CALIB()
//! Inital saturation of elements near river segments
{
	return rivEle_CALIB;
}

realtype setis_CALIB()
//! Calibration factor for maximum interception storage capacity
{
	return is_CALIB;
}
realtype setet0_CALIB()
//! Calibration factor cum Switch for ET0 (Evaporation from Cannopy)
{
	return et0_CALIB;
}
realtype setmf_CALIB()
//! Calibration factor for melt factor (rate)
{
	return mf_CALIB;
}
realtype settf_CALIB()
//! Calibration factor cum Switch for through-fall
{
	return tf_CALIB;
}


realtype setVic_CALIB()
//! Calibration fraction for Variable rate of infiltration
{
	return Vic_CALIB;
}
realtype setrivK_CALIB()
//! Calibration factor for hydraulic conductivity at river-element interface
{
	return rivK_CALIB;
}
realtype setKh_CALIB()
//! Calibration factor for horizontal hydraulic conductivity
{
	return Kh_CALIB;
}
realtype setRec_CALIB()
//! Calibration factor for recharge rate
{
	return Rec_CALIB;
}
realtype setet2_CALIB()
//! Calibration factor cum Switch for ET2
{
	return et2_CALIB;
}
realtype setet1_CALIB()
//! Calibration factor cum Switch for ET1
{
	return et1_CALIB;
}
realtype setsat_THRESH()
//! Saturation value above which macropore plays role
{
	return sat_THRESH;
}
realtype setmp_MULTFH()
//! Horizontal Hydraulic conductivity factor for macropore zones
{
	return mp_MULTFH;
}
realtype setmp_MULTFV()
//! Vertical Hydraulic conductivity factor for macropore zones
{
	return mp_MULTFV;
}
realtype setmpArea_CALIB()
//! Fraction of area governed by macropore behavior
{
	return mpArea_CALIB;
}
realtype setovl_THRESH_H(){
	return ovl_THRESH_H;
}
realtype setovl_THRESH_V(){
	return ovl_THRESH_V;
}
realtype setrzd_CALIB()
//! Uniform value for Root Zone Depth
{
	return rzd_CALIB;
}


realtype setroughEle_CALIB()
//! Calibration factor for the manning's roughness coefficient of elements
{
	return roughEle_CALIB;
}
realtype setroughRiv_CALIB()
//! Calibration factor for the manning's roughness coefficient of river segments
{
	return roughRiv_CALIB;
}
realtype setrivCoeff_CALIB()
//! Calibration factor for width of the river segments
{
	return rivCoeff_CALIB;
}
realtype setrivDepth_CALIB()
//! Calibration factor for depth of the river segments
{
	return rivDepth_CALIB;
}
realtype setalpha_CALIB()
//! Calibration factor for soil alpha paramter
{
	return alpha_CALIB;
}
realtype setset_MP()
//! Switch to turn Macropore behavior ON/OFF
{
	return set_MP;
}
realtype setlai_CALIB()
//! Calibration factor for Leaf Area Indices
{
	return lai_CALIB;
}
realtype setvegfrac_CALIB()
//! Calibration factor for vegetation fractions
{
	return vegfrac_CALIB;
}
realtype setalbedo_CALIB()
//! Calibration factor for albedo
{
	return albedo_CALIB;
}
