#ifndef PIHM_H
#define PIHM_H

 /******************************************************************************
 * File        : pihm.h                                                        *
 * Function    : define data structure and grobal variable                     *
 * Programmers : Yizhong Qu   @ Pennsylvania State Univeristy                  *
 *               Mukesh Kumar @ Pennsylvania State Univeristy                  *
 *               Gopal Bhatt  @ Pennsylvania State Univeristy                  *
 * Version     : v2.0 (July 10, 2007)                                          *
 *                                                                               *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 * This is an object oriented model. In this file, all the data structures and *
 * global variables are defined. To avoid using real globle variables, this    *
 * header file and pointers are passed.                                        *
 *                                                                             *
 * This code is free for users with research purpose only, if appropriate      *
 * citation is refered. However, there is no warranty in any format for this   *
 * product.                                                                    *
 *                                                                             *
 * For questions or comments, please contact the authors of the reference.     *
 * One who want to use it for other consideration may also contact Dr.Duffy    *
 * at cxd11@psu.edu.                                                           *
 *******************************************************************************/

//! @file pihm.h All the data structures and global variables are defined in this file

/* INCLUDE FILES */
#include <stdio.h>



/* SUNDIALS LIBRARY INCLUDE */
#include "sundials_types.h"
#include "nvector_serial.h"



/* Define Global Type */
//! Variable to store current time
float Tsteps;                 /* Variable to store current time */



/*******************************************************************************/
/*                       PIHM data structures                                  */
/*******************************************************************************/

/* Data Structure of an Triangular Element */
typedef struct element_type
//! Data Structure of an Triangular Element
{
    int index;                /**< Serial number of an element                   */
    int node[3];              /**< nodes of an element :: counter-clock-wise     */
    int nabr[3];              /**< neighbor i shares edge i (0: on boundary)     */

    realtype edge[3];         /**< edge[i] is formed by node[i] and node[i+1]    */
    realtype area;            /**< Area of an element                            */

    realtype x;               /**< x of centroid                                 */
    realtype y;               /**< y of centroid                                 */
    realtype zmin;            /**< z_min of centroid                             */
    realtype zmax;            /**< z_max of centroid                             */
    realtype NodeZmin;        /**< Z_min of any of the three nodes               */
    realtype NodeZmax;        /**< Z_max of any of the three nodes               */
    realtype NodeDist;        /**< Distance between Z_max and Z_min nodes        */

    realtype Ksat;            /**< Saturated Hydraulic Conductivity              */
    realtype Porosity;        /**< Porosity of the element                       */
    realtype Alpha;           /**<                                               */
    realtype Beta;            /**< Exponent                                      */
    realtype Sf;              /**< Friction Slope of the element                 */
    realtype RzD;             /**< Root Zone Depth                               */


    realtype LAImax;          /**< Maximum Leaf Area Index                       */
    realtype VegFrac;         /**< Vegitation Fraction                           */
    realtype Albedo;          /**< Albedo                                        */
    realtype Rs_ref;          /**<                                               */
    realtype Rmin;            /**< Stomatal Resistence                           */
    realtype Rough;           /**< Roughness (Manning's) Coefficient             */
    realtype windH;           /**< Height at which Wind Velocity is measured     */

    /* Class/Type associated with an Element */
    int soil;                 /**< soil type                                     */
    int LC;                   /**< Land Cover type                               */
    int IC;                   /**< initial condition type                        */
    int BC;                   /**< boundary type. 0:natural bc (no flow);
                                  +ve:Dirichlet BC; -ve:Nuemann BC               */
    int prep;                 /**< precipatation (forcing) type                  */
    int temp;                 /**< tempurature (forcing) type                    */
    int humidity;             /**< humidity type                                 */
    int WindVel;              /**< wind velocity type                            */
    int Rn;                   /**< net radiation input                           */
    int G;                    /**< radiation into ground                         */
    int pressure;             /**< pressure type                                 */
    int source;               /**< source (well) type                            */

} element;



/* Data Structure of an Node/Point */
typedef struct nodes_type
//! Data Structure of an Node/Point
{
    int index;                /**< Node Number                                   */

    realtype x;               /**< x coordinate of the node                      */
    realtype y;               /**< y coordinate of the node                      */
    realtype zmin;            /**< bedrock elevation of the node                 */
    realtype zmax;            /**< surface elevation  of the node                */

} nodes;


/* Data Structure of Initial Condition State of an Element */
typedef struct element_IC_type
//! Data Structure of Initial Condition State of an Element
{
    int index;

    realtype interception;    /**< Interception IC State (head) of the Element   */
    realtype snow;            /**< Snow IC State (head) of the Element           */
    realtype surf;            /**< Overland IC State (head) of the Element       */
    realtype unsat;           /**< Unsaturated IC State (head) of the Element    */
    realtype sat;             /**< Saturated IC State (head) of the Element      */

} element_IC;



/* Data Structure of a Soil Type */
typedef struct soils_type
//! Data Structure of a Soil Type
{
    int index;               /**< soil type class number                         */

    realtype Ksat;           /**< saturated soil conductivity                    */
    realtype SitaS;          /**< soil porosity                                  */
    realtype SitaR;          /**< soil moisture residual                         */
    realtype Alpha;          /**< soil curve parameter 1                         */
    realtype Beta;           /**< soil curve parameter 2                         */
    realtype Sf;             /**< surface slope of friction                      */
    realtype RzD;            /**< rootZone Depth                                 */

    int Macropore;           /**< 1: macropore; 0: regular soil                  */
    realtype base;           /**< base value                                     */
    realtype gama;           /**< amplifier factor                               */

    int Inf;                 /**< index of infiltration capacity type            */

} soils;



/* Data Structure of a Land Cover Type */
typedef struct lc_type
//! Data Structure of a Land Cover Type
{
    int index;               /**< land cover type index                          */

    realtype LAImax;         /**< max LAI                                        */
    realtype VegFrac;        /**< Canopy Fraction                                */
    realtype Albedo;         /**< Albedo                                         */
    realtype Rs_ref;         /**< reference stomaral resistance                  */
    realtype Rmin;           /**< Minimum stomatal resistance                    */
    realtype Rough;          /**< roughness (manning's) coefficient              */
} LC;



/* Data Structure of a River Segment */
typedef struct river_segment_type
//! Data Structure of a River Segment
{
    int index;                /**< Segment Number of the River                     */

    /* Topology of a River Segment */
    realtype x;               /**< x-cordinate of center of the river segment      */
    realtype y;               /**< y-cordinate of center of the river segment      */
    realtype zmin;            /**< bed elevation                                   */
    realtype zmax;            /**< bank elevation                                  */
    realtype depth;           /**< max depth of river segment                      */
    realtype Length;          /**< Length of the river segment                     */
    realtype Rough;

    /* Topological Relations of a River Segment */
    int FromNode;             /**< upstream Node Number                            */
    int ToNode;               /**< downstream Node Number                          */
    int down;                 /**< downstream segment Number                       */
    int LeftEle;              /**< Element at the Left of the segment              */
    int RightEle;             /**< Element at the Right of the segment             */
    int shape;                /**< shape type of the segment                       */
    int material;             /**< material type                                   */
    int IC;                   /**< IC type                                         */
    int BC;                   /**< BC type                                         */
    int reservoir;

} river_segment;



/* Data Structure of a River Shape */
typedef struct river_shape_type
//! Data Structure of a River Shape
{
	//!River Shape Type Number
    int index;                /**< River Shape Type Number                         */

    realtype width;           /**< assume rectangular shape                        */
    realtype depth;           /**< depth                                           */
    realtype bed;             /**< bed elevation                                   */
    int interpOrd;            /**< Interpolation order for river shape:
              Order= 1 (rectangle), 2(triangle), 3(quadratic) and 4(cubic)         */
    realtype coeff;           /**< Coefficient c in D = c*B/2                      */

} river_shape;



/* Data Structure of a River Material Type */
typedef struct river_material_type
//! Data Structure of a River Material Type
{
    int index;                /**< River Material Type Number                     */

    realtype Rough;           /**< Roughness (Mannin's) Coefficient               */
    realtype Sf;              /**< Friction Slope                                 */
    realtype Cwr;             /**< Discharge Coefficient of the edge              */

} river_material;



/* Data Structure of Initial Condition of the River Segment */
typedef struct river_IC_type
//! Data Structure of Initial Condition of the River Segment
{
    int index;                /**< IC State Type                                  */
    realtype value;           /**< Initial State (Head) of the River Segment      */

} river_IC;



/* Data Structure of Time Series Data Types */
typedef struct TSD_type
//! Data Structure of Time Series Data Type
{
    char name[50];            /**< Name of the Time Series :: Dummy Variable      */
    int index;                /**< Time Series Number                             */
    int length;               /**< Length of the Time Series                      */
    int iCounter;             /**< Current Time Series Access Pointer Index       */
    realtype **TS;            /**< 2D time series data type                       */

} TSD;



/* Model Data Structure */
typedef struct model_data_structure
//! Model (PIHM) Data Structure
{
    /* PIHM Mode Identifiers: Used to switch between several modes */
    int UnsatMode;               /**< Unsat Mode Identifier                       */
    int SurfMode;                /**< Surface Overland Mode Identifier            */
    int RivMode;                 /**< River Routing Mode Identifier               */

    /* Number of different model representation components */
    int NumEle;                  /**< Number of Elements in the model domain      */
    int NumNode;                 /**< Number of Nodes in the model domain         */
    int NumRiv;                  /**< Number of River Segs in the model domain    */

      /* Number of Unique Forcing TimeSeries within the model domain */
    int NumPrep;                 /**< Number of Precipitation TS                  */
    int NumTemp;                 /**< Number of Temperature TS                    */
    int NumHumidity;             /**< Number of Humidity TS                       */
    int NumWindVel;              /**< Number of Wind Velocity TS                  */
    int NumRn;                   /**< Number of Net Radiation TS                  */
    int NumG;                    /**< Number of Ground Heat TS                    */
    int NumP;                    /**< Number of Pressure TS                       */
    int NumSource;               /**< Number of Sources/Sinks TS                  */
    int NumMeltF;                /**< Number of Melt Factor : NumMeltF = 1        */

    /* Number of Spatial Attributes  of the Elements */
    int NumSoil;                 /**< Number of Soil Type                         */
    int NumRes;                  /**< Number of Reservoir/Dams                    */
    int NumInc;                  /**< Number of Infiltration Capacity             */
    int NumLC;                   /**< Number of Land Cover Index Data             */

    /* Number of Initial and Boundary condition TimeSeries for Elements */
    int Num1BC;                  /**< Number of Dirichlet BC                      */
    int Num2BC;                  /**< Number of Numann BC                         */
    int NumEleIC;                /**< Number of Element Initial Condtion          */

    /* Number of Spatial Attributes of River Segments */
    int NumRivShape;             /**< Number of River Shape                       */
    int NumRivMaterial;          /**< Number of River Bank Material               */

    /* Number of Initial and Boundary condition TimeSeries for Rivers */
    int NumRivIC;                /**< Number of River Initial Condition           */
    int NumRivBC;                /**< Number of River Boundary Condition          */

    /* Objects for feature (element/river) data types */
    element *Ele;                /**< Element (Triangle) Object Information       */
    river_segment *Riv;          /**< River Segment Information                   */

    /* Attributes of Element (triangle) objects in the model domain */
    nodes *Node;                 /**< Node Information                            */
    element_IC *Ele_IC;          /**< Element Initial Condtion                    */
    soils *Soil;                 /**< Soil Information                            */
    LC *LandC;                        /* Land Cover Information */

    /* Attributes of River (linear) objects in the model domain */
    river_shape *Riv_Shape;      /**< River Shape Information                     */
    river_material *Riv_Mat;     /**< River Bank Material Information             */
    river_IC *Riv_IC;            /**< River Initial Condition                     */

    /* Time Series Data in the model domain */
    TSD *TSD_Inc;                /**< Infiltration Capacity                       */
    TSD *TSD_LAI;                /**< Leaf Area Index Time Series Data            */
    TSD *TSD_DH;                 /**< Zero plane Displacement Height              */
    realtype *SIFactor;          /**< SIFactor :to calculate SIMax from LAI       */

    TSD *TSD_MeltF;              /**< Melt Factor for Temperature Index model     */

    TSD *TSD_EleBC;              /**< Element Boundary Condition                  */
    TSD *TSD_Prep;               /**< Precipitation Time Series Data              */
    TSD *TSD_Temp;               /**< Temperature Time Series Data                */
    TSD *TSD_Humidity;           /**< Humidity Time Series Data                   */
    TSD *TSD_WindVel;            /**< Wind Velocity Time Series Data              */
    TSD *TSD_Rn;                 /**< Net Radiation Time Series Data              */
    TSD *TSD_G;                  /**< Radiation into Ground Time Series Data      */
    TSD *TSD_Pressure;           /**< Pressure Time Series data                   */
    TSD *TSD_Source;             /**< Source (well) Time Series data              */

    realtype *WindH;             /**< Height at which wind velocity is observed   */

    TSD *TSD_Riv;                /**< River Related Time Series Data              */

    /* Storage for fluxes at Time = t */
    realtype **FluxSurf;         /**< Overland Flux between two elements          */
    realtype **FluxSub;          /**< Subsurface Flux between two elements        */
    realtype **FluxRiv;          /**< Flux between River Segs and Elements        */

    realtype *ElePrep;           /**< Rate of Prepicitation                       */
    realtype *Ele2IS;            /**< Rate of Interception                        */
    realtype *EleNetPrep;        /**< Net Rate of Precipitation                   */
    realtype *EleVic;            /**< Infiltration Rate                           */
    realtype *Recharge;          /**< Rate of Recharge (to GroundWater/Sat)       */
    realtype *EleSnow;           /**< Snow Accumulation/Storage                   */
    realtype *EleIS;             /**< Interception Storage                        */
    realtype *EleISmax;          /**< Maximum Interception Storage Capacity       */
    realtype *EleTF;             /**< Rate of Through Fall                        */
    realtype **EleET;            /**< Rate of Evapo-Transpiration                 */
    realtype Q;

} *Model_Data;



/* Control Data Structure */
typedef struct control_data_structure
//! Control Data Structure
{
    /* Model Control Option */
    int Verbose;                 /**< 1: Yea 0: Nay                               */
    int Debug;                   /**< 1: Yea 0: Nay                               */
    int int_type;                /**< Input Mode : 0: relax 1: typing 2:file      */
    int res_out;                 /**< 1: Yea 0: Nay                               */
    int flux_out;                /**< 1: Yea 0: Nay                               */
    int q_out;                   /**< 1: Yea 0: Nay                               */
    int etis_out;                /**< 1: Yea 0: Nay                               */

    /* Solver Control Options */
    int Solver;                  /**< Solver Type::                               */
    realtype abstol;             /**< Absolute Tolerance                          */
    realtype reltol;             /**< Relative Tolerance                          */
    realtype InitStep;           /**< Initial step size                           */
    realtype MaxStep;            /**< Maximum absolute step size                  */
    realtype ETStep;             /**< Absolute step size for ET Computation       */

    int GSType, MaxK;
    realtype delt;

    realtype StartTime;          /**< Simulation Start (Real) Time                */
    realtype EndTime;            /**< Simulation End (Real) Time                  */

    /* CVode is called at t = b * a^i */
    int outtype;
    realtype a;                  /**< Output Step Size Factor                     */
    realtype b;                  /**< Base Step Size                              */
    int NumSteps;                /**< Number of Step to be taken by CVode         */
    realtype *Tout;              /**< Array of Time at which State is computed    */

} Control_Data;



#endif
