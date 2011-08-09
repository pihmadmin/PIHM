/*******************************************************************************
 * File        : print.h                                                       *
 * Function    : defines identifiers for print.c file                          *
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

//! @file print.h define output file mode and which output paramater to print. Also, corresponding function declaration

/*    File Print Mode identifiers    */
#define TXT            1        /**< TXT is read as 1 */
#define CDF            2        /**< CDF is read as 2 */

/*    Fprint Control identifiers    */
#define YEA            1        /**< YEA is read as 1 */
#define NAY            0        /**< NAY is read as 0 */


/************************/
/* 1=> .txt Files       */
/* 2=> .netcdf Files    */
/************************/
#define FPRINT_MODE    CDF		/**< Specify output file mode: 1=.txt; 2=.nc */


//////////////////////////////////
#define ISState        YEA      /**< Output interception storage state? YEA:NAY */
#define SatState       YEA      /**< Output staturated zone state? YEA:NAY      */
#define UsatState      YEA      /**< Output unsaturated zone state? YEA:NAY     */
#define SurfState      YEA      /**< Output surface state? YEA:NAY              */

#define ISStateT       60       /**< Output mean interception storage at _ minute intervel   */
#define SatStateT      60       /**< Output mean saturated zone state at _ minute intervel   */
#define UsatStateT     60       /**< Output mean unsaturated zone state at _ minute intervel */
#define SurfStateT     60       /**< Output mean surface state at _ minute intervel          */
//////////////////////////////////
#define ET0            YEA      /**< Output evaporation rate from cannopy? YEA:NAY */
#define ET1            YEA      /**< Output ET1? YEA:NAY                           */
#define ET2            YEA      /**< Output ET2? YEA:NAY                           */
#define NetPpt         YEA      /**< Output net precipitation rate? YEA:NAY        */
#define Infil          YEA      /**< Output infiltration rate? YEA:NAY             */
#define RECHARGE       YEA      /**< Output recharge rate to ground water? YEA:NAY */

#define ET0T           60       /**< Output mean evaporation rate from cannopy at _ minute intervel   */
#define ET1T           60       /**< Output mean et1 at _ minute intervel                             */
#define ET2T           60       /**< Output mean et2 at _ minute intervel                             */
#define NetPptT        60       /**< Output mean net precipitation rate at _ minute intervel          */
#define InfilT         60       /**< Output infiltration rate at _ minute intervel                    */
#define RECHARGET      60       /**< Output recharge rate to groundwater at _ minute intervel         */
//////////////////////////////////
#define RivHead        YEA      /**< Output head of river segments? YEA:NAY */

#define RivHeadT       60       /**< Output river head at _ minute intervel         */
//////////////////////////////////
#define RivFlow        YEA      /**< Output outflow from river segments? YEA:NAY   */
#define RivBase        YEA      /**< Output baseflow to river segments? YEA:NAY    */
#define RivSurf        YEA      /**< Output surfaceflow to river segments? YEA:NAY */

#define RivFlowT       60       /**< Output outflow from river segments at _ minute intervel   */
#define RivBaseT       60       /**< Output baseflow to river segments at _ minute intervel    */
#define RivSurfT       60       /**< Output surfaceflow to river segments at _ minute intervel */
//////////////////////////////////
#define FluxSurf       YEA
#define FluxSat        YEA

#define FluxSurfT      60
#define FluxSatT       60
//////////////////////////////////


/* Function Prototypes */
void FPrintInit(Model_Data);                                        /* Initialize variables for FPrint           */

void printIS(Model_Data, FILE *, realtype);                         /* Print Interception Storage in TXT mode    */
void printIScdf(Model_Data, int, int, realtype);                    /* Print Interception Storage in CDF mode    */

void printSatState(Model_Data, N_Vector, FILE *, realtype);         /* Print Saturated State in TXT mode         */
void printSatStatecdf(Model_Data, N_Vector, int, int, realtype);    /* Print Saturated State in CDF mode         */
void printUsatState(Model_Data, N_Vector, FILE *, realtype);        /* Print Unsaturated State in TXT mode       */
void printUsatStatecdf(Model_Data, N_Vector, int, int, realtype);   /* Print Unsaturated State in CDF mode       */
void printSurfState(Model_Data, N_Vector, FILE *, realtype);        /* Print Overland State in TXT mode          */
void printSurfStatecdf(Model_Data, N_Vector, int, int, realtype);   /* Print Overland State in CDF mode          */

void printET0(Model_Data, FILE *, realtype);                        /* Print Evaporation from Veg in TXT mode    */
void printET0cdf(Model_Data, int, int, realtype);                   /* Print Evaporation from Veg in CDF mode    */
void printET1(Model_Data, FILE *, realtype);                        /* Print Evaporation from Sur in TXT mode    */
void printET1cdf(Model_Data, int, int, realtype);                   /* Print Evaporation from Sur in CDF mode    */
void printET2(Model_Data, FILE *, realtype);                        /* Print Evaporation from Grd in TXT mode    */
void printET2cdf(Model_Data, int, int, realtype);                   /* Print Evaporation from Grd in CDF mode    */
void printNetPpt(Model_Data, FILE *, realtype);                     /* Print Net Precipitation in TXT mode       */
void printNetPptcdf(Model_Data, int, int, realtype);                /* Print Net Precipitation in CDF mode       */
void printInfil(Model_Data, FILE *, realtype);                      /* Print Variable Infiltration in TXT mode   */
void printInfilcdf(Model_Data, int, int, realtype);                 /* Print Variable Infiltration in CDF mode   */
void printRecharge(Model_Data, FILE *, realtype);                   /* Print Recharge to GW in TXT mode          */
void printRechargecdf(Model_Data, int, int, realtype);              /* Print Recharge to GW in CDF mode          */

void printRiverFlow(Model_Data, N_Vector, FILE *, realtype);        /* Print outflow from river segin TXT mode   */
void printRiverFlowcdf(Model_Data, N_Vector, int, int, realtype);   /* Print outflow from river segin CDF mode   */
void printRiverBase(Model_Data, FILE *, realtype);                  /* Print Base flow to river seg in TXT mode  */
void printRiverBasecdf(Model_Data, int, int, realtype);             /* Print Base flow to river seg in CDF mode  */
void printRiverSurf(Model_Data, FILE *, realtype);                  /* Print over flow to river seg in TXT mode  */
void printRiverSurfcdf(Model_Data, int, int, realtype);             /* Print over flow to river seg in CDF mode  */
void printRiverHead(Model_Data, N_Vector, FILE *, realtype);        /* Print River State in TXT mode             */
void printRiverHeadcdf(Model_Data, N_Vector, int, int, realtype);   /* Print River State in CDF mode             */

void FPrintCloseAll(void);
