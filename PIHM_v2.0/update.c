/*******************************************************************************
 * File        : update.c                                                      *
 * Function    : Update variables during simulation 		   	       *
 * Version     : Nov, 2007 (2.0)                                               *
 * Developer of PIHM2.0:        Mukesh Kumar (muk139@psu.edu)                  *
 * Developer of PIHM1.0:        Yizhong Qu                                     *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0............................*
 * a) New File for resetting TS index counter                                  *
 * b) Can be used in refinement/(de)refinement simulation or for any other     *
 *    temporal update of model parameters				       *  
 * Acknowledgement: Thanks to Bhatt, G. for idenfication of inefficiency in    *
 * Interpolation funcn.							       *
 *******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "nvector_serial.h"
#include "sundialstypes.h"   
#include "pihm.h"      

#define EPSILON 0.05
#define UNIT_C 1440     /* 60*24 for calculation of yDot in m/min units while forcing is in m/day. */

realtype Interpolation(TSD *Data, realtype t);

void update(realtype t, void *DS)
	{
  	int k;
  
  	Model_Data MD;
  
  	MD = (Model_Data)DS;
        for(k=0; k<MD->NumPrep; k++)
        	{
                while(MD->TSD_Prep[k].iCounter < MD->TSD_Prep[k].length && t/(UNIT_C) > MD->TSD_Prep[k].TS[MD->TSD_Prep[k].iCounter+1][0])
			{
                        MD->TSD_Prep[k].iCounter++;
                       	}
       		}
       	for(k=0; k<MD->NumTemp; k++)
             	{
              	while(MD->TSD_Temp[k].iCounter < MD->TSD_Temp[k].length && t/(UNIT_C) > MD->TSD_Temp[k].TS[MD->TSD_Temp[k].iCounter+1][0])
               		{
                 	MD->TSD_Temp[k].iCounter++;
               		}
               	}
      	for(k=0; k<MD->NumHumidity; k++)
             	{
            	while(MD->TSD_Humidity[k].iCounter < MD->TSD_Humidity[k].length && t/(UNIT_C) > MD->TSD_Humidity[k].TS[MD->TSD_Humidity[k].iCounter+1][0])
              		{
            		MD->TSD_Humidity[k].iCounter++;
              		}
          	}
     	for(k=0; k<MD->NumWindVel; k++)
        	{   
         	while(MD->TSD_WindVel[k].iCounter < MD->TSD_WindVel[k].length && t/(UNIT_C) > MD->TSD_WindVel[k].TS[MD->TSD_WindVel[k].iCounter+1][0])
        		{
        		MD->TSD_WindVel[k].iCounter++;
        		}
      		}
 
    	for(k=0; k<MD->NumRn; k++)
     		{                               
 		while(MD->TSD_Rn[k].iCounter < MD->TSD_Rn[k].length && t/(UNIT_C) > MD->TSD_Rn[k].TS[MD->TSD_Rn[k].iCounter+1][0])
     			{
      			MD->TSD_Rn[k].iCounter++;
    			}
   		}
   	for(k=0; k<MD->NumG; k++)
     		{  
      		while(MD->TSD_G[k].iCounter < MD->TSD_G[k].length && t/(UNIT_C) > MD->TSD_G[k].TS[MD->TSD_G[k].iCounter+1][0])
          		{
          		MD->TSD_G[k].iCounter++;
           		}
    		}
 	for(k=0; k<MD->NumP; k++)  
   		{     
   		while(MD->TSD_Pressure[k].iCounter < MD->TSD_Pressure[k].length && t/(UNIT_C) > MD->TSD_Pressure[k].TS[MD->TSD_Pressure[k].iCounter+1][0])
        		{
       			MD->TSD_Pressure[k].iCounter++;
     			}
    		}
 	for(k=0; k<MD->NumLC; k++)  
    		{    
    		while(MD->TSD_LAI[k].iCounter < MD->TSD_LAI[k].length && t/(UNIT_C) > MD->TSD_LAI[k].TS[MD->TSD_LAI[k].iCounter+1][0])
      			{
      			MD->TSD_LAI[k].iCounter++;
      			}
   		}
  	for(k=0; k<MD->NumLC; k++)
    		{     
   		while(MD->TSD_RL[k].iCounter < MD->TSD_RL[k].length && t/(UNIT_C) > MD->TSD_RL[k].TS[MD->TSD_RL[k].iCounter+1][0])
     			{
      			MD->TSD_RL[k].iCounter++;
      			}
   		}
 	for(k=0; k<MD->NumMeltF; k++)  
     		{           
  		while(MD->TSD_MeltF[k].iCounter < MD->TSD_MeltF[k].length && t/(UNIT_C) > MD->TSD_MeltF[k].TS[MD->TSD_MeltF[k].iCounter+1][0])
      			{
      			MD->TSD_MeltF[k].iCounter++;
      			}
 		}
 	for(k=0; k<MD->NumSource; k++)  
    		{     
   		while(MD->TSD_Source[k].iCounter < MD->TSD_Source[k].length && t/(UNIT_C) > MD->TSD_Source[k].TS[MD->TSD_Source[k].iCounter+1][0])
  			{
   			MD->TSD_Source[k].iCounter++;
  			}
		}

	}
