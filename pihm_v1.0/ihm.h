/*******************************************************************************
 * File        : ihm.h                                                         *
 * Function    : define data structure and grobal variable                     *
 * Programmer  : Yizhong Qu @ Pennsylvania State Univeristy                    *
 * Version     : May, 2004 (1.0)                                               * 
 *-----------------------------------------------------------------------------*
 *                                                                             *
 * This is a object oriented model. In this file, all the data structures and  *
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

#include "sundialstypes.h"
#include "nvector_serial.h"

/* Define Global Type */

typedef struct element_type 
{
  int index;  
  int node[3];        /* anti-clock-wise */
  int nabr[3];        /* neighbor i shares edge i (0: on boundary) */
    
  realtype edge[3];   /* edge i is from node i to node i+1 */
  realtype area;      /* area of element */
  
  realtype x;         /* x of centroid */
  realtype y;         /* y of centroid */
  realtype zmin;      /* z_min of centroid */
  realtype zmax;      /* z_max of centroid */
  
  realtype Ksat;
  realtype Porosity;
  realtype Alpha;
  realtype Beta;
  realtype Sf;
  realtype Rough;
    
  int soil;           /* soil type */
  int LAI;            /* LAI type  */
  int IC;             /* initial condition type */
  int BC;             /* boundary type. 0:natural bc (no flow); 1:Dirichlet BC; 2:Nuemann BC */
  int prep;           /* precipatation (forcing) type */
  int temp;           /* tempurature (forcing) type   */
  int humidity;       /* humidity type */
  int WindVel;        /* wind velocity type  */
  int Rn;             /* net radiation input */
  int G;              /* radiation into ground */
  int pressure;       /* pressure type */
  int source;         /* source (well) type */
  
} element;

typedef struct nodes_type
{
  int index; 
     
  realtype x;          /* x coordinate */
  realtype y;          /* y coordinate */
  realtype zmin;       /* z bed rock elevation */
  realtype zmax;       /* z surface elevation  */
  
} nodes;

typedef struct element_IC_type
{
  int index;
  
  realtype interception;
  realtype surf;
  realtype unsat;
  realtype sat;
  
} element_IC;

typedef struct soils_type
{
  int index;           /* index */
  
  realtype Ksat;       /* saturated soil conductivity */
  realtype SitaS;      /* soil porosity */
  realtype SitaR;      /* soil moisture residual */
  realtype Alpha;      /* soil curve parameter 1 */
  realtype Beta;       /* soil curve parameter 2 */
  realtype Sf;         /* surface slope of friction */
  realtype Rough;      /* surface roughness factor  */
  
  int Macropore;       /* 1: macropore; 0: regular soil */
  realtype base;       /* base value */
  realtype gama;       /* amplifier factor */
  
  int Inf;             /* index of infiltration capacity type */
  
} soils;

typedef struct river_segment_type
{
  int index;
  
  realtype x;          /* x of river segment */
  realtype y;  
  realtype zmin;       /* bed elevation  */
  realtype zmax;       /* bank elevation */
  realtype depth;      /* max depth */
  realtype Length;
  realtype Rough;
  
  int FromNode;
  int ToNode;
  int down;            /* down stream segment */
  int LeftEle;
  int RightEle;
  int shape;           /* shape type    */
  int material;        /* material type */
  int IC;              /* IC type */
  int BC;              /* BC type */
  int reservoir;
  
} river_segment;

typedef struct river_shape_type
{
  int index;           
  realtype width;      /* assume rectangular shape */
  realtype depth;      /* depth */
  realtype bed;        /* bed elevation */
  
} river_shape;

typedef struct river_material_type
{
  int index;
  realtype Rough;
  realtype Sf;
  realtype Cwr;
  
} river_material;

typedef struct river_IC_type
{
  int index;
  realtype value;
  
} river_IC;

typedef struct TSD_type
{
  char name[5];
  int index;
  int length;
  realtype **TS;     /* 2D time series data */
  
} TSD;

typedef struct model_data_structure          /* Model_data definition */
{
  int UnsatMode;               /* Unsat Mode */
  int SurfMode;                /* Surface Overland Mode */
  int RivMode;                 /* River Routing Mode */
    
  int NumEle;                  /* Number of Elements */
  int NumNode;                 /* Number of Nodes    */
  int NumRiv;                  /* Number of Rivere Segments */
  
  int NumPrep;                 /* Number of Precipatation   */
  int NumTemp;                 /* Number of Temperature     */
  int NumHumidity;             /* Number of Humidity        */
  int NumWindVel;              /* Number of Wind Velocity   */
  int NumRn;                   /* Number of Net Radiation   */
  int NumG;                    /* Number of Ground Heat     */
  int NumP;                    /* Number of Pressure        */
  int NumSource;               /* Number of Source          */
  
  int NumSoil;                 /* Number of Soils           */
  int NumRes;                  /* Number of Reservoir       */
  int NumInc;                  /* Number of Infiltration Capacity  */
  int NumLAI;                  /* Number of Leaves Area Index Data */
  
  int Num1BC;                  /* Number of Dirichlet BC    */
  int Num2BC;                  /* Number of Numann BC       */
  int NumEleIC;                /* Number of Element Initial Condtion */
  
  int NumRivShape;             /* Number of River Shape     */
  int NumRivMaterial;          /* Number of River Bank Material */
  int NumRivIC;                /* Number of River Initial Condition  */
  int NumRivBC;                /* Number of River Boundary Condition */
  
  element *Ele;                /* Store Element Information  */
  nodes *Node;                 /* Store Node Information     */
  element_IC *Ele_IC;          /* Store Element Initial Condtion */
  soils *Soil;                 /* Store Soil Information     */
  
  river_segment *Riv;          /* Store River Segment Information */
  river_shape *Riv_Shape;      /* Store River Shape Information   */
  river_material *Riv_Mat;     /* Store River Bank Material Information */
  river_IC *Riv_IC;            /* Store River Initial Condition   */

  TSD *TSD_Inc;                /* Infiltration Capacity Time Series Data */
  TSD *TSD_LAI;                /* Leaves Area Index Time Series Data     */
  realtype *SIFactor;          /* SIFactor is used to calculate SIMax from LAI */

  TSD *TSD_EleBC;              /* Element Boundary Condition Time Series Data  */
  TSD *TSD_Riv;                /* River Related Time Series Data  */
  TSD *TSD_Prep;               /* RainFall Time Series Data       */
  TSD *TSD_Temp;               /* Temperature Time Series Data    */
  TSD *TSD_Humidity;           /* Humidity Time Series Data       */
  TSD *TSD_WindVel;            /* Wind Velocity Time Series Data  */
  TSD *TSD_Rn;                 /* Net Radiation Time Series Data  */
  TSD *TSD_G;                  /* Radiation into Ground Time Series Data */
  TSD *TSD_Pressure;           /* Pressure Time Series data       */
  TSD *TSD_Source;             /* Source (well) Time Series data  */

  realtype **FluxSurf;     /* Overland Flux   */
  realtype **FluxSub;      /* Subsurface Flux */
  realtype **FluxRiv;      /* River Segement Flux */

  realtype *ElePrep;
  realtype *Ele2IS;
  realtype *EleNetPrep;
  realtype *EleVic;
  realtype *Recharge;
  realtype *EleIS;
  realtype *EleISmax;
  realtype *EleETP;
  realtype **EleET;
  realtype Q; 
  
} *Model_Data;

typedef struct control_data_structure
{
  int Verbose;
  int Debug;
  
  int Solver;
  int NumSteps;
  
  int res_out;
  int flux_out;
  int q_out;
  int etis_out;
  
  int int_type;
  
  realtype abstol;
  realtype reltol;
  realtype InitStep;
  realtype MaxStep;
  realtype ETStep;
  
  int GSType, MaxK;
  realtype delt;
  
  realtype StartTime;
  realtype EndTime;
  
  int outtype;
  realtype a;
  realtype b;
  
  realtype *Tout;  
  
} Control_Data;


void PrintModelData(Model_Data);
void PrintEle(Model_Data);
void PrintEleAtt(Model_Data);
void PrintNode(Model_Data);
void PrintSoil(Model_Data);
void PrintTS(TSD *, int);
void PrintRiv(Model_Data);
void PrintForcing(Model_Data);
void PrintDY(Model_Data, N_Vector, N_Vector);
void PrintY(Model_Data, N_Vector, realtype);
void FPrintYheader(FILE *, Model_Data);
void FPrintY(Model_Data, N_Vector, realtype, FILE *);
void FPrintFlux(Model_Data, realtype, FILE *);
void FPrintETISheader(FILE *, Model_Data);
void FPrintETIS(Model_Data, realtype, FILE *);
void FPrintQ(Model_Data, realtype, FILE *);
void PrintVerbose(int, realtype, long int iopt[], realtype ropt[]);
void PrintFarewell(Control_Data, long int iopt[], realtype ropt[], realtype, realtype);
void PrintFinalStats(long int iopt[], realtype ropt[]);
void FPrintFarewell(Control_Data, FILE *, long int iopt[], realtype ropt[], realtype, realtype);
void FPrintFinalStats(FILE *, long int iopt[], realtype ropt[]);

