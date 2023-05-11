/*
 *   ex_sfun_ctrl_VSI_AFE.c Simple C-MEX S-function for function call.
 *
 *   ABSTRACT:
 *     The purpose of this sfunction is to call a simple legacy
 *     function during simulation:
 *
 *        void epwm1_isr(double y1[1], double y2[1], double y3[1], double y4[1], double y5[1], double y6[1], double y7[1], double y8[1], double y9[1], double u1, double u2, double u3, double u4, double u5, double u6, double u7, double u8, double u9, double u10, double u11, double u12, double u13, double u14, double u15, double u16, double u17, double u18, double u19, double u20, double u21)
 *
 *   Simulink version           : 8.1 (R2013a) 13-Feb-2013
 *   C source code generated on : 11-Aug-2014 14:59:01

 * THIS S-FUNCTION IS GENERATED BY THE LEGACY CODE TOOL AND MAY NOT WORK IF MODIFIED

 */

/*
   %%%-MATLAB_Construction_Commands_Start
   def = legacy_code('initialize');
   def.SFunctionName = 'ex_sfun_ctrl_VSI_AFE';
   def.OutputFcnSpec = 'void epwm1_isr(double y1[1], double y2[1], double y3[1], double y4[1], double y5[1], double y6[1], double y7[1], double y8[1], double y9[1], double u1, double u2, double u3, double u4, double u5, double u6, double u7, double u8, double u9, double u10, double u11, double u12, double u13, double u14, double u15, double u16, double u17, double u18, double u19, double u20, double u21)';
   def.HeaderFiles = {'math.h'};
   def.SourceFiles = {'epwm1_isr.c'};
   def.SampleTime = 5e-05;
   legacy_code('sfcn_cmex_generate', def);
   legacy_code('compile', def);
   %%%-MATLAB_Construction_Commands_End
 */

/*
 * Must specify the S_FUNCTION_NAME as the name of the S-function.
 */
#define S_FUNCTION_NAME                ex_sfun_ctrl_VSI_AFE
#define S_FUNCTION_LEVEL               2

/*
 * Need to include simstruc.h for the definition of the SimStruct and
 * its associated macro definitions.
 */
#include "simstruc.h"

/*
 * Specific header file(s) required by the legacy code function.
 */
#include "math.h"

/* Function: mdlInitializeSizes ===========================================
 * Abstract:
 *    The sizes information is used by Simulink to determine the S-function
 *    block's characteristics (number of inputs, outputs, states, etc.).
 */
static void mdlInitializeSizes(SimStruct *S)
{
  /* Number of expected parameters */
  ssSetNumSFcnParams(S, 0);

  /*
   * Set the number of pworks.
   */
  ssSetNumPWork(S, 0);

  /*
   * Set the number of dworks.
   */
  if (!ssSetNumDWork(S, 0))
    return;

  /*
   * Set the number of input ports.
   */
  if (!ssSetNumInputPorts(S, 21))
    return;

  /*
   * Configure the input port 1
   */
  ssSetInputPortDataType(S, 0, SS_DOUBLE);
  ssSetInputPortWidth(S, 0, 1);
  ssSetInputPortComplexSignal(S, 0, COMPLEX_NO);
  ssSetInputPortDirectFeedThrough(S, 0, 1);
  ssSetInputPortAcceptExprInRTW(S, 0, 1);
  ssSetInputPortOverWritable(S, 0, 1);
  ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
  ssSetInputPortRequiredContiguous(S, 0, 1);

  /*
   * Configure the input port 2
   */
  ssSetInputPortDataType(S, 1, SS_DOUBLE);
  ssSetInputPortWidth(S, 1, 1);
  ssSetInputPortComplexSignal(S, 1, COMPLEX_NO);
  ssSetInputPortDirectFeedThrough(S, 1, 1);
  ssSetInputPortAcceptExprInRTW(S, 1, 1);
  ssSetInputPortOverWritable(S, 1, 1);
  ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
  ssSetInputPortRequiredContiguous(S, 1, 1);

  /*
   * Configure the input port 3
   */
  ssSetInputPortDataType(S, 2, SS_DOUBLE);
  ssSetInputPortWidth(S, 2, 1);
  ssSetInputPortComplexSignal(S, 2, COMPLEX_NO);
  ssSetInputPortDirectFeedThrough(S, 2, 1);
  ssSetInputPortAcceptExprInRTW(S, 2, 1);
  ssSetInputPortOverWritable(S, 2, 1);
  ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
  ssSetInputPortRequiredContiguous(S, 2, 1);

  /*
   * Configure the input port 4
   */
  ssSetInputPortDataType(S, 3, SS_DOUBLE);
  ssSetInputPortWidth(S, 3, 1);
  ssSetInputPortComplexSignal(S, 3, COMPLEX_NO);
  ssSetInputPortDirectFeedThrough(S, 3, 1);
  ssSetInputPortAcceptExprInRTW(S, 3, 1);
  ssSetInputPortOverWritable(S, 3, 1);
  ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
  ssSetInputPortRequiredContiguous(S, 3, 1);

  /*
   * Configure the input port 5
   */
  ssSetInputPortDataType(S, 4, SS_DOUBLE);
  ssSetInputPortWidth(S, 4, 1);
  ssSetInputPortComplexSignal(S, 4, COMPLEX_NO);
  ssSetInputPortDirectFeedThrough(S, 4, 1);
  ssSetInputPortAcceptExprInRTW(S, 4, 1);
  ssSetInputPortOverWritable(S, 4, 1);
  ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
  ssSetInputPortRequiredContiguous(S, 4, 1);

  /*
   * Configure the input port 6
   */
  ssSetInputPortDataType(S, 5, SS_DOUBLE);
  ssSetInputPortWidth(S, 5, 1);
  ssSetInputPortComplexSignal(S, 5, COMPLEX_NO);
  ssSetInputPortDirectFeedThrough(S, 5, 1);
  ssSetInputPortAcceptExprInRTW(S, 5, 1);
  ssSetInputPortOverWritable(S, 5, 1);
  ssSetInputPortOptimOpts(S, 5, SS_REUSABLE_AND_LOCAL);
  ssSetInputPortRequiredContiguous(S, 5, 1);

  /*
   * Configure the input port 7
   */
  ssSetInputPortDataType(S, 6, SS_DOUBLE);
  ssSetInputPortWidth(S, 6, 1);
  ssSetInputPortComplexSignal(S, 6, COMPLEX_NO);
  ssSetInputPortDirectFeedThrough(S, 6, 1);
  ssSetInputPortAcceptExprInRTW(S, 6, 1);
  ssSetInputPortOverWritable(S, 6, 1);
  ssSetInputPortOptimOpts(S, 6, SS_REUSABLE_AND_LOCAL);
  ssSetInputPortRequiredContiguous(S, 6, 1);

  /*
   * Configure the input port 8
   */
  ssSetInputPortDataType(S, 7, SS_DOUBLE);
  ssSetInputPortWidth(S, 7, 1);
  ssSetInputPortComplexSignal(S, 7, COMPLEX_NO);
  ssSetInputPortDirectFeedThrough(S, 7, 1);
  ssSetInputPortAcceptExprInRTW(S, 7, 1);
  ssSetInputPortOverWritable(S, 7, 1);
  ssSetInputPortOptimOpts(S, 7, SS_REUSABLE_AND_LOCAL);
  ssSetInputPortRequiredContiguous(S, 7, 1);

  /*
   * Configure the input port 9
   */
  ssSetInputPortDataType(S, 8, SS_DOUBLE);
  ssSetInputPortWidth(S, 8, 1);
  ssSetInputPortComplexSignal(S, 8, COMPLEX_NO);
  ssSetInputPortDirectFeedThrough(S, 8, 1);
  ssSetInputPortAcceptExprInRTW(S, 8, 1);
  ssSetInputPortOverWritable(S, 8, 1);
  ssSetInputPortOptimOpts(S, 8, SS_REUSABLE_AND_LOCAL);
  ssSetInputPortRequiredContiguous(S, 8, 1);

  /*
   * Configure the input port 10
   */
  ssSetInputPortDataType(S, 9, SS_DOUBLE);
  ssSetInputPortWidth(S, 9, 1);
  ssSetInputPortComplexSignal(S, 9, COMPLEX_NO);
  ssSetInputPortDirectFeedThrough(S, 9, 1);
  ssSetInputPortAcceptExprInRTW(S, 9, 1);
  ssSetInputPortOverWritable(S, 9, 1);
  ssSetInputPortOptimOpts(S, 9, SS_REUSABLE_AND_LOCAL);
  ssSetInputPortRequiredContiguous(S, 9, 1);

  /*
   * Configure the input port 11
   */
  ssSetInputPortDataType(S, 10, SS_DOUBLE);
  ssSetInputPortWidth(S, 10, 1);
  ssSetInputPortComplexSignal(S, 10, COMPLEX_NO);
  ssSetInputPortDirectFeedThrough(S, 10, 1);
  ssSetInputPortAcceptExprInRTW(S, 10, 1);
  ssSetInputPortOverWritable(S, 10, 1);
  ssSetInputPortOptimOpts(S, 10, SS_REUSABLE_AND_LOCAL);
  ssSetInputPortRequiredContiguous(S, 10, 1);

  /*
   * Configure the input port 12
   */
  ssSetInputPortDataType(S, 11, SS_DOUBLE);
  ssSetInputPortWidth(S, 11, 1);
  ssSetInputPortComplexSignal(S, 11, COMPLEX_NO);
  ssSetInputPortDirectFeedThrough(S, 11, 1);
  ssSetInputPortAcceptExprInRTW(S, 11, 1);
  ssSetInputPortOverWritable(S, 11, 1);
  ssSetInputPortOptimOpts(S, 11, SS_REUSABLE_AND_LOCAL);
  ssSetInputPortRequiredContiguous(S, 11, 1);

  /*
   * Configure the input port 13
   */
  ssSetInputPortDataType(S, 12, SS_DOUBLE);
  ssSetInputPortWidth(S, 12, 1);
  ssSetInputPortComplexSignal(S, 12, COMPLEX_NO);
  ssSetInputPortDirectFeedThrough(S, 12, 1);
  ssSetInputPortAcceptExprInRTW(S, 12, 1);
  ssSetInputPortOverWritable(S, 12, 1);
  ssSetInputPortOptimOpts(S, 12, SS_REUSABLE_AND_LOCAL);
  ssSetInputPortRequiredContiguous(S, 12, 1);

  /*
   * Configure the input port 14
   */
  ssSetInputPortDataType(S, 13, SS_DOUBLE);
  ssSetInputPortWidth(S, 13, 1);
  ssSetInputPortComplexSignal(S, 13, COMPLEX_NO);
  ssSetInputPortDirectFeedThrough(S, 13, 1);
  ssSetInputPortAcceptExprInRTW(S, 13, 1);
  ssSetInputPortOverWritable(S, 13, 1);
  ssSetInputPortOptimOpts(S, 13, SS_REUSABLE_AND_LOCAL);
  ssSetInputPortRequiredContiguous(S, 13, 1);

  /*
   * Configure the input port 15
   */
  ssSetInputPortDataType(S, 14, SS_DOUBLE);
  ssSetInputPortWidth(S, 14, 1);
  ssSetInputPortComplexSignal(S, 14, COMPLEX_NO);
  ssSetInputPortDirectFeedThrough(S, 14, 1);
  ssSetInputPortAcceptExprInRTW(S, 14, 1);
  ssSetInputPortOverWritable(S, 14, 1);
  ssSetInputPortOptimOpts(S, 14, SS_REUSABLE_AND_LOCAL);
  ssSetInputPortRequiredContiguous(S, 14, 1);

  /*
   * Configure the input port 16
   */
  ssSetInputPortDataType(S, 15, SS_DOUBLE);
  ssSetInputPortWidth(S, 15, 1);
  ssSetInputPortComplexSignal(S, 15, COMPLEX_NO);
  ssSetInputPortDirectFeedThrough(S, 15, 1);
  ssSetInputPortAcceptExprInRTW(S, 15, 1);
  ssSetInputPortOverWritable(S, 15, 1);
  ssSetInputPortOptimOpts(S, 15, SS_REUSABLE_AND_LOCAL);
  ssSetInputPortRequiredContiguous(S, 15, 1);

  /*
   * Configure the input port 17
   */
  ssSetInputPortDataType(S, 16, SS_DOUBLE);
  ssSetInputPortWidth(S, 16, 1);
  ssSetInputPortComplexSignal(S, 16, COMPLEX_NO);
  ssSetInputPortDirectFeedThrough(S, 16, 1);
  ssSetInputPortAcceptExprInRTW(S, 16, 1);
  ssSetInputPortOverWritable(S, 16, 1);
  ssSetInputPortOptimOpts(S, 16, SS_REUSABLE_AND_LOCAL);
  ssSetInputPortRequiredContiguous(S, 16, 1);

  /*
   * Configure the input port 18
   */
  ssSetInputPortDataType(S, 17, SS_DOUBLE);
  ssSetInputPortWidth(S, 17, 1);
  ssSetInputPortComplexSignal(S, 17, COMPLEX_NO);
  ssSetInputPortDirectFeedThrough(S, 17, 1);
  ssSetInputPortAcceptExprInRTW(S, 17, 1);
  ssSetInputPortOverWritable(S, 17, 1);
  ssSetInputPortOptimOpts(S, 17, SS_REUSABLE_AND_LOCAL);
  ssSetInputPortRequiredContiguous(S, 17, 1);

  /*
   * Configure the input port 19
   */
  ssSetInputPortDataType(S, 18, SS_DOUBLE);
  ssSetInputPortWidth(S, 18, 1);
  ssSetInputPortComplexSignal(S, 18, COMPLEX_NO);
  ssSetInputPortDirectFeedThrough(S, 18, 1);
  ssSetInputPortAcceptExprInRTW(S, 18, 1);
  ssSetInputPortOverWritable(S, 18, 1);
  ssSetInputPortOptimOpts(S, 18, SS_REUSABLE_AND_LOCAL);
  ssSetInputPortRequiredContiguous(S, 18, 1);

  /*
   * Configure the input port 20
   */
  ssSetInputPortDataType(S, 19, SS_DOUBLE);
  ssSetInputPortWidth(S, 19, 1);
  ssSetInputPortComplexSignal(S, 19, COMPLEX_NO);
  ssSetInputPortDirectFeedThrough(S, 19, 1);
  ssSetInputPortAcceptExprInRTW(S, 19, 1);
  ssSetInputPortOverWritable(S, 19, 1);
  ssSetInputPortOptimOpts(S, 19, SS_REUSABLE_AND_LOCAL);
  ssSetInputPortRequiredContiguous(S, 19, 1);

  /*
   * Configure the input port 21
   */
  ssSetInputPortDataType(S, 20, SS_DOUBLE);
  ssSetInputPortWidth(S, 20, 1);
  ssSetInputPortComplexSignal(S, 20, COMPLEX_NO);
  ssSetInputPortDirectFeedThrough(S, 20, 1);
  ssSetInputPortAcceptExprInRTW(S, 20, 1);
  ssSetInputPortOverWritable(S, 20, 1);
  ssSetInputPortOptimOpts(S, 20, SS_REUSABLE_AND_LOCAL);
  ssSetInputPortRequiredContiguous(S, 20, 1);

  /*
   * Set the number of output ports.
   */
  if (!ssSetNumOutputPorts(S, 9))
    return;

  /*
   * Configure the output port 1
   */
  ssSetOutputPortDataType(S, 0, SS_DOUBLE);
  ssSetOutputPortWidth(S, 0, 1);
  ssSetOutputPortComplexSignal(S, 0, COMPLEX_NO);
  ssSetOutputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
  ssSetOutputPortOutputExprInRTW(S, 0, 0);

  /*
   * Configure the output port 2
   */
  ssSetOutputPortDataType(S, 1, SS_DOUBLE);
  ssSetOutputPortWidth(S, 1, 1);
  ssSetOutputPortComplexSignal(S, 1, COMPLEX_NO);
  ssSetOutputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
  ssSetOutputPortOutputExprInRTW(S, 1, 0);

  /*
   * Configure the output port 3
   */
  ssSetOutputPortDataType(S, 2, SS_DOUBLE);
  ssSetOutputPortWidth(S, 2, 1);
  ssSetOutputPortComplexSignal(S, 2, COMPLEX_NO);
  ssSetOutputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
  ssSetOutputPortOutputExprInRTW(S, 2, 0);

  /*
   * Configure the output port 4
   */
  ssSetOutputPortDataType(S, 3, SS_DOUBLE);
  ssSetOutputPortWidth(S, 3, 1);
  ssSetOutputPortComplexSignal(S, 3, COMPLEX_NO);
  ssSetOutputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
  ssSetOutputPortOutputExprInRTW(S, 3, 0);

  /*
   * Configure the output port 5
   */
  ssSetOutputPortDataType(S, 4, SS_DOUBLE);
  ssSetOutputPortWidth(S, 4, 1);
  ssSetOutputPortComplexSignal(S, 4, COMPLEX_NO);
  ssSetOutputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
  ssSetOutputPortOutputExprInRTW(S, 4, 0);

  /*
   * Configure the output port 6
   */
  ssSetOutputPortDataType(S, 5, SS_DOUBLE);
  ssSetOutputPortWidth(S, 5, 1);
  ssSetOutputPortComplexSignal(S, 5, COMPLEX_NO);
  ssSetOutputPortOptimOpts(S, 5, SS_REUSABLE_AND_LOCAL);
  ssSetOutputPortOutputExprInRTW(S, 5, 0);

  /*
   * Configure the output port 7
   */
  ssSetOutputPortDataType(S, 6, SS_DOUBLE);
  ssSetOutputPortWidth(S, 6, 1);
  ssSetOutputPortComplexSignal(S, 6, COMPLEX_NO);
  ssSetOutputPortOptimOpts(S, 6, SS_REUSABLE_AND_LOCAL);
  ssSetOutputPortOutputExprInRTW(S, 6, 0);

  /*
   * Configure the output port 8
   */
  ssSetOutputPortDataType(S, 7, SS_DOUBLE);
  ssSetOutputPortWidth(S, 7, 1);
  ssSetOutputPortComplexSignal(S, 7, COMPLEX_NO);
  ssSetOutputPortOptimOpts(S, 7, SS_REUSABLE_AND_LOCAL);
  ssSetOutputPortOutputExprInRTW(S, 7, 0);

  /*
   * Configure the output port 9
   */
  ssSetOutputPortDataType(S, 8, SS_DOUBLE);
  ssSetOutputPortWidth(S, 8, 1);
  ssSetOutputPortComplexSignal(S, 8, COMPLEX_NO);
  ssSetOutputPortOptimOpts(S, 8, SS_REUSABLE_AND_LOCAL);
  ssSetOutputPortOutputExprInRTW(S, 8, 0);

  /*
   * Register reserved identifiers to avoid name conflict
   */
  if (ssRTWGenIsCodeGen(S) || ssGetSimMode(S)==SS_SIMMODE_EXTERNAL) {
    /*
     * Register reserved identifier for OutputFcnSpec
     */
    ssRegMdlInfo(S, "epwm1_isr", MDL_INFO_ID_RESERVED, 0, 0, ssGetPath(S));
  }

  /*
   * This S-function can be used in referenced model simulating in normal mode.
   */
  ssSetModelReferenceNormalModeSupport(S, MDL_START_AND_MDL_PROCESS_PARAMS_OK);

  /*
   * Set the number of sample time.
   */
  ssSetNumSampleTimes(S, 1);

  /*
   * All options have the form SS_OPTION_<name> and are documented in
   * matlabroot/simulink/include/simstruc.h. The options should be
   * bitwise or'd together as in
   *   ssSetOptions(S, (SS_OPTION_name1 | SS_OPTION_name2))
   */
  ssSetOptions(S,
               SS_OPTION_USE_TLC_WITH_ACCELERATOR |
               SS_OPTION_CAN_BE_CALLED_CONDITIONALLY |
               SS_OPTION_EXCEPTION_FREE_CODE |
               SS_OPTION_WORKS_WITH_CODE_REUSE |
               SS_OPTION_SFUNCTION_INLINED_FOR_RTW |
               SS_OPTION_DISALLOW_CONSTANT_SAMPLE_TIME);
}

/* Function: mdlInitializeSampleTimes =====================================
 * Abstract:
 *    This function is used to specify the sample time(s) for your
 *    S-function. You must register the same number of sample times as
 *    specified in ssSetNumSampleTimes.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
  ssSetSampleTime(S, 0, (real_T)5e-05);
  ssSetOffsetTime(S, 0, (real_T)0);

#if defined(ssSetModelReferenceSampleTimeDisallowInheritance)

  ssSetModelReferenceSampleTimeDisallowInheritance(S);

#endif

}

/* Function: mdlOutputs ===================================================
 * Abstract:
 *    In this function, you compute the outputs of your S-function
 *    block. Generally outputs are placed in the output vector(s),
 *    ssGetOutputPortSignal.
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
  /*
   * Get access to Parameter/Input/Output/DWork/size information
   */
  real_T *u1 = (real_T *) ssGetInputPortSignal(S, 0);
  real_T *u2 = (real_T *) ssGetInputPortSignal(S, 1);
  real_T *u3 = (real_T *) ssGetInputPortSignal(S, 2);
  real_T *u4 = (real_T *) ssGetInputPortSignal(S, 3);
  real_T *u5 = (real_T *) ssGetInputPortSignal(S, 4);
  real_T *u6 = (real_T *) ssGetInputPortSignal(S, 5);
  real_T *u7 = (real_T *) ssGetInputPortSignal(S, 6);
  real_T *u8 = (real_T *) ssGetInputPortSignal(S, 7);
  real_T *u9 = (real_T *) ssGetInputPortSignal(S, 8);
  real_T *u10 = (real_T *) ssGetInputPortSignal(S, 9);
  real_T *u11 = (real_T *) ssGetInputPortSignal(S, 10);
  real_T *u12 = (real_T *) ssGetInputPortSignal(S, 11);
  real_T *u13 = (real_T *) ssGetInputPortSignal(S, 12);
  real_T *u14 = (real_T *) ssGetInputPortSignal(S, 13);
  real_T *u15 = (real_T *) ssGetInputPortSignal(S, 14);
  real_T *u16 = (real_T *) ssGetInputPortSignal(S, 15);
  real_T *u17 = (real_T *) ssGetInputPortSignal(S, 16);
  real_T *u18 = (real_T *) ssGetInputPortSignal(S, 17);
  real_T *u19 = (real_T *) ssGetInputPortSignal(S, 18);
  real_T *u20 = (real_T *) ssGetInputPortSignal(S, 19);
  real_T *u21 = (real_T *) ssGetInputPortSignal(S, 20);
  real_T *y1 = (real_T *) ssGetOutputPortSignal(S, 0);
  real_T *y2 = (real_T *) ssGetOutputPortSignal(S, 1);
  real_T *y3 = (real_T *) ssGetOutputPortSignal(S, 2);
  real_T *y4 = (real_T *) ssGetOutputPortSignal(S, 3);
  real_T *y5 = (real_T *) ssGetOutputPortSignal(S, 4);
  real_T *y6 = (real_T *) ssGetOutputPortSignal(S, 5);
  real_T *y7 = (real_T *) ssGetOutputPortSignal(S, 6);
  real_T *y8 = (real_T *) ssGetOutputPortSignal(S, 7);
  real_T *y9 = (real_T *) ssGetOutputPortSignal(S, 8);

  /*
   * Call the legacy code function
   */
  epwm1_isr( y1, y2, y3, y4, y5, y6, y7, y8, y9, *u1, *u2, *u3, *u4, *u5, *u6,
            *u7, *u8, *u9, *u10, *u11, *u12, *u13, *u14, *u15, *u16, *u17, *u18,
            *u19, *u20, *u21);
}

/* Function: mdlTerminate =================================================
 * Abstract:
 *    In this function, you should perform any actions that are necessary
 *    at the termination of a simulation.
 */
static void mdlTerminate(SimStruct *S)
{
}

/*
 * Required S-function trailer
 */
#ifdef MATLAB_MEX_FILE
# include "simulink.c"
#else
# include "cg_sfun.h"
#endif
