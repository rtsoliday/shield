/* Copyright 2016 by Michael Borland and Argonne National Laboratory,
 * all rights reserved.
 */
/* program: shield.c 
 * purpose: A C port of the venerable SHIELD11 program
 *
 * Michael Borland, 2016
 */
#include "mdb.h"
#include "SDDS.h"
#include "scan.h"
#include "rpn.h"
#include "namelist.h"

#include "shield.h"

#define USAGE "shield <inputFile> -dataFile=<filename>\n\n\
inputFile      List of namelist commands for a particular shielding calculation.\n\
-dataFile      Give the name of the file that contains material properties used in calculations.\n\
-describeInput Describe namelist input.\n\n\
Program by M. Borland, July 2016. Based on SHIELD11 (SLAC).\n"

#define CLOPT_DATAFILE 0
#define CLOPT_DESCRIBE_INPUT 1
#define N_CLOPTS 2
char *clOption[N_CLOPTS] = {
  "datafile", "describeinput",
};

#define RUN_COMMAND 0
#define TARGET_COMMAND 1
#define PRIMARY_SHIELD_COMMAND 2
#define SECONDARY_SHIELD_COMMAND 3
#define N_COMMANDS 4

char *command[N_COMMANDS] = {
  "run", "target", "primary_shield", "secondary_shield",
};

typedef struct {
  long nMaterials;
  char **material;
  int32_t *Z;      /* atomic number */
  double *A;    /* atomic mass in g/mole */
  double *rho;  /* density in g/cm^3 */
  double *X0;   /* radiation length in g/cm^2 */
  double *R0;   /* Molier length in g/cm^2 */
  /* Mean free path in g/cm^2 */
  double *XmfpGRNs;  /* giant resonance neutrons */
  double *XmfpMIDs;  /* medium energy neutrons */
  double *XmfpHENs;  /* high-energy neutrons */
  double *XmfpGamD;      /* direct gamma */
  double *XmfpGamI;      /* indirect gamma */
} MATERIAL_DATA;

typedef struct {
  long materialID;
  double thickness_cm;
} SHIELD;


void processRunCommand(NAMELIST_TEXT *nltext);
void processTargetCommand(NAMELIST_TEXT *nltext, MATERIAL_DATA *mData);
void processPrimaryShieldCommand(NAMELIST_TEXT *nltext, SHIELD **shield, long *nShields, MATERIAL_DATA *materialData);
void processSecondaryShieldCommand(NAMELIST_TEXT *nltext, SHIELD **shield, long *nShields, MATERIAL_DATA *materialData);
int readMaterialDataFile(char *dataFile, MATERIAL_DATA *materialData);
void reportValidMaterials(char *message, char *material, MATERIAL_DATA *mData);
void setUpOutputFile(SDDS_DATASET *SDDSout, RUN_STRUCT *run, char *inputFile, char *mDataFile,
                     TARGET_STRUCT *target, PRIMARY_SHIELD_STRUCT *primary, SHIELD *shield, long nShields, MATERIAL_DATA *materialData);
void runShieldingCase(RUN_STRUCT *run, TARGET_STRUCT *target, SHIELD *shield, long nShields, double distance, double angleDeg, MATERIAL_DATA *mData, SDDS_DATASET *SDDSout);
double computeSourceTerm(char *component, TARGET_STRUCT *target, MATERIAL_DATA *mData, double energyGeV, double angle);

void link_date(void)
{
  fprintf(stdout, "Link date: %s %s, SVN revision: %s\n", __DATE__, __TIME__, SVN_VERSION);
  fflush(stdout);
}

int main(int argc, char **argv)
{
  FILE *fp_in;
  short target_seen, primary_shield_seen, secondary_shield_seen;
  char s[16834];
  NAMELIST_TEXT namelist_text;
  long i;
  char *inputFile, *dataFile;
  int i_arg;
  SCANNED_ARG *scanned;
  MATERIAL_DATA materialData;
  SHIELD *shield = NULL;
  long nShields = 0;
  SDDS_DATASET SDDSout;
  
  if ((argc = scanargsg(&scanned, argc, argv))<2) {
    link_date();
    bomb(NULL, USAGE);
  }
  
  inputFile = dataFile = NULL;
  
  for (i_arg=1; i_arg<argc; i_arg++) {
    if (scanned[i_arg].arg_type==OPTION) {
      /* process options here */
      switch (match_string(scanned[i_arg].list[0], clOption, N_CLOPTS, 0)) {
      case CLOPT_DATAFILE:
        if (scanned[i_arg].n_items!=2)
          bomb("invalid -dataFile syntax", NULL);
        dataFile = scanned[i_arg].list[1];
        break;
      case CLOPT_DESCRIBE_INPUT:
        show_namelists_fields(stdout, namelist_pointer, namelist_name, n_namelists);
        if (argc==2) exit(0);
        break;
      default:
        printf("option %s not recognized\n", scanned[i_arg].list[0]);
        bomb("unknown option given", USAGE);
        break;
      }
    }
    else {
      /* filename argument */
      if (!inputFile)
        inputFile = scanned[i_arg].list[0];
      else
        bomb("too many filenames", NULL);
    }
  }

  if (!(fp_in = fopen(inputFile, "r"))) 
    bomb("problem reading input file", NULL);
  
if (!dataFile || !strlen(dataFile))
    bomb("no material data file given", USAGE);
  
  if (!readMaterialDataFile(dataFile, &materialData))
    bomb("problem reading shield data file", NULL);
  
  target_seen = primary_shield_seen = secondary_shield_seen = 0;

  while (get_namelist(s, 16384, fp_in)) {
    scan_namelist(&namelist_text, s);
    switch (match_string(namelist_text.group_name, command, N_COMMANDS, EXACT_MATCH)) {
    case RUN_COMMAND:
      if (!target_seen)
        bomb("No target defined", NULL);
      if (!primary_shield_seen)
        bomb("No primary shield defined", NULL);
      processRunCommand(&namelist_text);
      setUpOutputFile(&SDDSout, &run_struct, inputFile, dataFile, &target_struct, &primary_shield_struct, shield, nShields, &materialData);
      runShieldingCase(&run_struct, &target_struct, shield, nShields, primary_shield_struct.distance_cm, primary_shield_struct.angle_deg, &materialData, &SDDSout);
      break;
    case TARGET_COMMAND:
      target_seen = 1;
      processTargetCommand(&namelist_text, &materialData);
      break;
    case PRIMARY_SHIELD_COMMAND:
      if (primary_shield_seen)
        bomb("Primary shield already seen", NULL);
      primary_shield_seen = 1;
      processPrimaryShieldCommand(&namelist_text, &shield, &nShields, &materialData);
      break;
    case SECONDARY_SHIELD_COMMAND:
      if (!primary_shield_seen)
        bomb("Please give the primary shield before any secondary shields", NULL);
      secondary_shield_seen = 1;
      processSecondaryShieldCommand(&namelist_text, &shield, &nShields, &materialData);
      break;
    default:
      printf("unknown namelist %s given.  Known namelists are:\n", namelist_text.group_name);
      for (i=0; i<N_COMMANDS; i++)
        printf("%s\n", command[i]);
      exit(1);
      break;
    }
  }
  return 0;
}

void setUpOutputFile(SDDS_DATASET *SDDSout, RUN_STRUCT *run, char *inputFile, char *mDataFile,
                     TARGET_STRUCT *target, PRIMARY_SHIELD_STRUCT *primary, SHIELD *shield, long nShields,
                     MATERIAL_DATA *mData)
{
  char buffer[16384];
  long i;
  
  sprintf(buffer, "shield output for input file %s with material data file %s", 
          inputFile, mDataFile);
  if (!SDDS_InitializeOutput(SDDSout, SDDS_BINARY, 1, buffer, NULL, run->output)) {
    SDDS_SetError("Unable to initialized output file");
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!SDDS_DefineSimpleParameter(SDDSout, "SVNVersion", NULL, SDDS_STRING) ||
      !SDDS_DefineSimpleParameter(SDDSout, "primaryShieldDistance", "cm", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(SDDSout, "primaryShieldAngle", "deg", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(SDDSout, "primaryShieldMaterial", NULL, SDDS_STRING) ||
      !SDDS_DefineSimpleParameter(SDDSout, "primaryShieldThickness", "cm", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(SDDSout, "targetMaterial", NULL, SDDS_STRING) ||
      !SDDS_DefineSimpleParameter(SDDSout, "targetLength", "cm", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleParameter(SDDSout, "targetAttenuation", NULL, SDDS_SHORT) ||
      !SDDS_DefineSimpleParameter(SDDSout, "targetRadius", "cm", SDDS_DOUBLE)) {
    SDDS_SetError("Unable to define parameters (1)");
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  for (i=1; i<nShields; i++) {
    sprintf(buffer, "secondaryShield%ldThickness", i);
    if (!SDDS_DefineSimpleParameter(SDDSout, buffer, "cm", SDDS_DOUBLE)) {
      SDDS_SetError("Unable to define parameters (2)");
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    sprintf(buffer, "secondaryShield%ldMaterial", i);
    if (!SDDS_DefineSimpleParameter(SDDSout, buffer, NULL, SDDS_STRING)) {
      SDDS_SetError("Unable to define parameters (3)");
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }
  if (!SDDS_DefineSimpleColumn(SDDSout, "theta", "deg", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "x", "cm", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "GRNDose", "rem/(kW*h)", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "MIDDose", "rem/(kW*h)", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "HENDose", "rem/(kW*h)", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "neutronDose", "rem/(kW*h)", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "GamIDose", "rem/(kW*h)", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "GamDDose", "rem/(kW*h)", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "gammaDose", "rem/(kW*h)", SDDS_DOUBLE) ||
      !SDDS_DefineSimpleColumn(SDDSout, "totalDose", "rem/(kW*h)", SDDS_DOUBLE)) {
    SDDS_SetError("Unable to define columns");
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if (!SDDS_WriteLayout(SDDSout) || !SDDS_StartPage(SDDSout, run->n_theta) ||
      !SDDS_SetParameters(SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                         "SVNVersion", SVN_VERSION,
                         "primaryShieldDistance", primary->distance_cm,
                         "primaryShieldAngle", primary->angle_deg,
                         "primaryShieldMaterial", primary->material,
                         "primaryShieldThickness", primary->thickness_cm,
                         "targetMaterial", target->material,
                         "targetLength", target->length_cm,
                         "targetAttenuation", (short)target->include_attenuation,
                         "targetRadius", target->radius_cm,
                         NULL)) {
    SDDS_SetError("Unable to set parameter values (1)");
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  for (i=1; i<nShields; i++) {
    sprintf(buffer, "secondaryShield%ldThickness", i);
    if (!SDDS_SetParameters(SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                            buffer, shield[i].thickness_cm, NULL)) {
      SDDS_SetError("Unable to set parameter values (2)");
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    sprintf(buffer, "secondaryShield%ldMaterial", i);
    if (!SDDS_SetParameters(SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE,
                            buffer, mData->material[shield[i].materialID], NULL)) {
      SDDS_SetError("Unable to set parameter values (2)");
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
  }
}


#define COMP_GRN 0
#define COMP_MID 1
#define COMP_HEN 2
#define COMP_GAMD 3
#define COMP_GAMI 4
#define N_COMPONENTS 5
static char *componentName[N_COMPONENTS] = {
  "GRN", "MID", "HEN", "GamD", "GamI",
};

void runShieldingCase(RUN_STRUCT *run, TARGET_STRUCT *target, 
                      SHIELD *shield, long nShields,  /* includes primary and secondary shields */
                      double distance, double angleDeg, /* distance and angle of shielding */
                      MATERIAL_DATA *mData, SDDS_DATASET *SDDSout)
{
  long iTheta, iComp;
  double theta, dTheta;
  double source[N_COMPONENTS], dose[N_COMPONENTS];
  double amtDeg, cosAMT;
  char buffer[16384];
  
  /* Scan theta */
  if (run->n_theta>1)
    dTheta = (run->theta_max_deg - run->theta_min_deg)/(run->n_theta-1);
  else
    dTheta = 0;
  printf("%7s %11s %11s %11s %11s %11s %11s %11s %11s %11s %11s\n",
         "theta", 
         "GRNsrc", "MIDsrc", "HENsrc", "GAMDsrc", "GAMIsrc",
         "GRN", "MID", "HEN", "GAMD", "GAMI");
  for (iTheta=0; iTheta<run->n_theta; iTheta++) {
    double neutronDose, gammaDose;
    theta = run->theta_min_deg + dTheta*iTheta;
    amtDeg = angleDeg - theta;
    cosAMT = cos(amtDeg*PI/180);
    if (amtDeg>90 || amtDeg<-90 || cosAMT==0) {
      printf("Warning: skipping shield angle (%le) and detector angle (%le), which differ by more than 90 deg\n", angleDeg, theta);
      continue;
    }
    printf("%7.3lf ", theta);
    for (iComp=0; iComp<N_COMPONENTS; iComp++) {
      source[iComp] 
        = computeSourceTerm(componentName[iComp], target, mData, run->beam_energy_GeV, theta);
      printf("%11.3e ", source[iComp]);
    }
    for (iComp=0; iComp<N_COMPONENTS; iComp++) {
      long sh;
      dose[iComp] = source[iComp]/sqr((distance+shield[0].thickness_cm)/cosAMT);   /* in rem/electron at outside of primary shield */
      for (sh=0; sh<nShields; sh++) {
        long m;
        m = shield[sh].materialID;
        switch (iComp) {
        case COMP_GRN:
          dose[iComp] *= exp(-shield[sh].thickness_cm/cosAMT*mData->rho[m]/mData->XmfpGRNs[m]);
          break;
        case COMP_MID:
          dose[iComp] *= exp(-shield[sh].thickness_cm/cosAMT*mData->rho[m]/mData->XmfpMIDs[m]);
          break;
        case COMP_HEN:
          dose[iComp] *= exp(-shield[sh].thickness_cm/cosAMT*mData->rho[m]/mData->XmfpHENs[m]);
          break;
        case COMP_GAMD:
          dose[iComp] *= exp(-shield[sh].thickness_cm/cosAMT*mData->rho[m]/mData->XmfpGamD[m]);
          break;
        case COMP_GAMI:
          dose[iComp] *= exp(-shield[sh].thickness_cm/cosAMT*mData->rho[m]/mData->XmfpGamI[m]);
          break;
        default:
          bomb("unknown component", NULL);
          break;
        }
      }
      printf("%11.3e ", dose[iComp]);
      sprintf(buffer, "%sDose", componentName[iComp]);
      if (!SDDS_SetRowValues(SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, iTheta,
                             buffer, dose[iComp]*2.25e16/run->beam_energy_GeV, NULL)) {
        SDDS_SetError("Unable to set dose column");
        SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
      }
    }
    neutronDose = gammaDose = 0;
    for (iComp=0; iComp<N_COMPONENTS; iComp++) {
      switch (iComp) {
      case COMP_GRN:
      case COMP_MID:
      case COMP_HEN:
        neutronDose += dose[iComp];
        break;
      default:
        gammaDose += dose[iComp];
        break;
      }
    }        
    if (!SDDS_SetRowValues(SDDSout, SDDS_SET_BY_NAME|SDDS_PASS_BY_VALUE, iTheta,
                           "theta", theta, "x", (distance+shield[0].thickness_cm)*tan(amtDeg*PI/180),
                           "neutronDose", neutronDose*2.25e16/run->beam_energy_GeV,
                           "gammaDose", gammaDose*2.25e16/run->beam_energy_GeV,
                           "totalDose", (neutronDose+gammaDose)*2.25e16/run->beam_energy_GeV,
                           NULL)) {
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
    }
    printf("\n");
  }
  if (!SDDS_WritePage(SDDSout) || !SDDS_Terminate(SDDSout)) {
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
}

double computeSourceTerm(char *component, TARGET_STRUCT *target, MATERIAL_DATA *mData, double energyGeV, double angle)
{
  long m, iComp;
  double lengthGPerSC, lengthX0, radiusGPerSC, radiusR0, radRelax;
  double SltSor, Sorc, cosTheta, sinTheta;
  double E_HEN[23] = {
              0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21,
              0.22, 0.23, 0.24, 0.25, 0.27, 0.30, 0.35,
              0.40, 0.45, 0.50, 0.55, 0.60, 0.70, 0.80,
              0.90, 1.00
              };
  double CS_HEN[23] = {
           0.00444,0.00711, 0.0111, 0.0156, 0.0222, 0.0298, 0.0382,
            0.0489, 0.0547, 0.0622, 0.0711, 0.0889,  0.116,  0.162,
             0.211,  0.276,  0.338,  0.404,  0.502,  0.601,  0.711,
             0.813,   1.00
           } ;

  double lengthNeutron = 17.332196; /* Minimum target length for neutrons (r.l.) */
  double radiusNeutron = 3.736411;  /* Minimum target radius for neutrons (Mol.units) */
  double thresholdHEN = 0.15 ;      /* Threshold energy (GeV) for HEN production */
  double lengthGamma =  0.01;       /* Minimum target length for photons (r.l.) */
  double radiusGamma = 1.189869;    /* Minimum target radius for photons */
  
  if ((m=match_string(target->material, mData->material, mData->nMaterials, EXACT_MATCH))<0) {
    reportValidMaterials("Invalid material for target", target_struct.material, mData);
    bomb("Invalid material used for target", NULL);
  }
  if ((iComp=match_string(component, componentName, N_COMPONENTS, EXACT_MATCH))<0) 
    bomb("invalid component in computeSourceTerm. Should never happen!", NULL);

  lengthGPerSC = target->length_cm*mData->rho[m]; /* length in gm/cm^2 */
  lengthX0 = lengthGPerSC/mData->X0[m];
  radiusGPerSC = target->radius_cm*mData->rho[m]; /* radius in gm/cm^2 */
  radiusR0 = radiusGPerSC/mData->R0[m];
  radRelax = radiusGPerSC/mData->XmfpGamD[m];

  cosTheta = cos(angle*PI/180);
  sinTheta = sin(angle*PI/180);
  SltSor = 0;
  
  if (target->include_attenuation && iComp!=COMP_GAMD) {
    if (lengthX0<lengthNeutron) {
      printf("Error for %s: Target length (%g in X0 units) is less than minimum value of %g\n",
             componentName[iComp], lengthX0, lengthNeutron);
      exit(1);
    }
    if (radiusR0<radiusNeutron) {
      printf("Error: Target radius (%g in X0 units) is less than minimum value of %g\n",
              radiusR0, radiusNeutron);
      exit(1);
    }
    if (fabs(angle)<=90) {
      double CritTar, CritCor, SltTar, SltCor;
      /* Find critical angle and SltTar for cylindrical target */
      CritTar = atan(target->radius_cm/target->length_cm)*180/PI;
      if (fabs(angle)>CritTar) 
        SltTar = radiusGPerSC/sinTheta;
      else
        SltTar = lengthGPerSC/cosTheta;
      
      /* Find critical angle and SltCor for cylindrical core */
      CritCor = atan(radiusNeutron*mData->R0[m]/(lengthNeutron*mData->X0[m]))*180/PI;
      if (fabs(angle)>=CritCor) 
        SltCor = radiusNeutron*mData->R0[m]/sinTheta;
      else
        SltCor = lengthNeutron*mData->X0[m]/cosTheta;

      SltSor = SltTar - SltCor;
    }
  }
  
  switch (iComp) {
  case COMP_GRN:
    Sorc = 4.93*pow(mData->Z[m], 0.662)*(energyGeV*1e-11);
    if (target->include_attenuation)
      Sorc *= exp(-SltSor/mData->XmfpGRNs[m]);
    break;
  case COMP_MID:
    Sorc = 43.9*energyGeV*1e-11*pow(mData->A[m], -0.37)/(1-0.75*cosTheta);
    if (energyGeV<=0.5) {
      Sorc *= 1.6*pow(energyGeV, 1.5);
    } else if (energyGeV>0.5 && energyGeV<1) {
      Sorc *= (0.566 + 0.434*(energyGeV-0.5)/0.5);
    }
    if (target->include_attenuation)
      Sorc *= exp(-SltSor/mData->XmfpMIDs[m]);
    break;
  case COMP_HEN:
  case COMP_GAMI:
    if (energyGeV<=thresholdHEN)
      Sorc = 0;
    else
      Sorc = 13.7*pow(mData->A[m], -0.65)/ipow(1-0.72*cosTheta, 2)*energyGeV*1e-11;
    if (energyGeV<1.0) {
      unsigned long interpCode;
      Sorc *= interpolate(CS_HEN, E_HEN, 23, energyGeV, NULL, NULL, 1, &interpCode, 1);
    }
    if (target->include_attenuation)
      Sorc *= exp(-SltSor/mData->XmfpHENs[m]);
    if (iComp==COMP_GAMI)
      Sorc *= 0.267;
    break;
  case COMP_GAMD:
    if (lengthX0 <= lengthGamma) {
      printf("Error: Target is too short for direct photon model\n");
      exit(1);
    }
    if (radRelax < radiusGamma) {
      printf("Error: Target radius is too small for direct photon model\n");
      exit(1);
    }
    Sorc = 1e6*energyGeV*exp(-lengthGPerSC/mData->XmfpGamD[m])*exp(-0.959*sqrt(fabs(angle)));
    if (fabs(angle)<90) {
      Sorc += 683*exp(-radRelax)*exp(-fabs(angle)/72.2);
    } else {
      Sorc += 683*exp(-radiusGamma)*exp(-fabs(angle)/72.2);
    }
    Sorc *= energyGeV*1e-11;
    break;
  default:
    printf("invalid component code %ld seen in computeSourceTerm---seek professional help!\n", iComp);
    exit(1);
    break;
  }
  return Sorc;
}



void processRunCommand(NAMELIST_TEXT *nltext)
{
  process_namelist(&run, nltext);
  print_namelist(stdout, &run);
  if (!run_struct.output || !strlen(run_struct.output))
    bomb("no output file named", NULL);
  if (run_struct.beam_energy_GeV<0.05)
    bomb("beam energy too low", NULL);
  if (run_struct.theta_min_deg>run_struct.theta_max_deg)
    bomb("theta_min_deg > theta_max_deg", NULL);
  if (run_struct.theta_min_deg<0)
    bomb("theta_min_deg<0", NULL);
  if (run_struct.theta_max_deg>180)
    bomb("theta_max_deg>180", NULL);
  if (run_struct.n_theta<1)
    bomb("n_theta<1", NULL);
}

void processTargetCommand(NAMELIST_TEXT *nltext, MATERIAL_DATA *mData)
{
  long mCode;
  process_namelist(&target, nltext);
  print_namelist(stdout, &target);
  if (!target_struct.material)
    bomb("No material named for target", NULL);
  if ((mCode=match_string(target_struct.material, mData->material, mData->nMaterials, EXACT_MATCH))<0) {
    reportValidMaterials("Invalid material for target", target_struct.material, mData);
    bomb("Invalid material used for target", NULL);
  }
  /* To do: check length and radius against radiation and Moliere lengths, respectively */
  if (target_struct.length_cm<=0)
    bomb("Target length must be positive", NULL);
  if (target_struct.radius_cm<=0)
    bomb("Target radius must be positive", NULL);
}

void reportValidMaterials(char *message, char *material, MATERIAL_DATA *mData)
{
  long i;
  printf("%s: %s\n", message, material);
  printf("Valid materials are:\n");
  for (i=0; i<mData->nMaterials; i++)
    printf("    \"%s\"\n", mData->material[i]);
}

void processPrimaryShieldCommand(NAMELIST_TEXT *nltext, SHIELD **shield, long *nShields, MATERIAL_DATA *mData)
{
  long mCode;
  process_namelist(&primary_shield, nltext);
  print_namelist(stdout, &primary_shield);

  if (!primary_shield_struct.material)
    bomb("No material named for primary shield", NULL);
  if ((mCode=match_string(primary_shield_struct.material, mData->material, mData->nMaterials, EXACT_MATCH))<0) {
    reportValidMaterials("Invalid material for primary shield", primary_shield_struct.material, mData);
    bomb("Invalid material used for primary shield", NULL);
  }
  if (primary_shield_struct.thickness_cm<0)
    bomb("Primary shield thickness must be non-negative", NULL);
  *shield = SDDS_Realloc(*shield, sizeof(**shield)*(*nShields+1));
  (*shield)[*nShields].materialID = mCode;
  (*shield)[*nShields].thickness_cm = primary_shield_struct.thickness_cm;
  *nShields += 1;
}

void processSecondaryShieldCommand(NAMELIST_TEXT *nltext, SHIELD **shield, long *nShields, MATERIAL_DATA *mData)
{
  long mCode;
  process_namelist(&secondary_shield, nltext);
  print_namelist(stdout, &secondary_shield);

  if (!secondary_shield_struct.material)
    bomb("No material named for secondary shield", NULL);
  if ((mCode=match_string(secondary_shield_struct.material, mData->material, mData->nMaterials, EXACT_MATCH))<0) {
    reportValidMaterials("Invalid material for secondary shield", secondary_shield_struct.material, mData);
    bomb("Invalid material used for secondary shield", NULL);
  }
  if (secondary_shield_struct.thickness_cm<=0)
    bomb("Secondary shield thickness must be positive", NULL);
  *shield = SDDS_Realloc(*shield, sizeof(**shield)*(*nShields+1));
  (*shield)[*nShields].materialID = mCode;
  (*shield)[*nShields].thickness_cm = secondary_shield_struct.thickness_cm;
  *nShields += 1;
}

int readMaterialDataFile(char *dataFile, MATERIAL_DATA *materialData)
{
  SDDS_DATASET SDDSin;
  long i;
  
  if (!SDDS_InitializeInput(&SDDSin, dataFile) || !SDDS_ReadPage(&SDDSin)) {
      SDDS_SetError("Unable to read material data file");
      SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }
  if ((materialData->nMaterials=SDDS_RowCount(&SDDSin))<=0) 
    bomb("No data in material data file", NULL);

  if (!(materialData->material = SDDS_GetColumn(&SDDSin, "Material")) ||
      !(materialData->Z = SDDS_GetColumnInLong(&SDDSin, "Z")) ||
      !(materialData->A = SDDS_GetColumnInDoubles(&SDDSin, "A")) ||
      !(materialData->rho = SDDS_GetColumnInDoubles(&SDDSin, "rho")) ||
      !(materialData->X0 = SDDS_GetColumnInDoubles(&SDDSin, "X0")) ||
      !(materialData->R0 = SDDS_GetColumnInDoubles(&SDDSin, "R0")) ||
      !(materialData->XmfpGRNs = SDDS_GetColumnInDoubles(&SDDSin, "XmfpGRNs")) ||
      !(materialData->XmfpMIDs = SDDS_GetColumnInDoubles(&SDDSin, "XmfpMIDs")) ||
      !(materialData->XmfpHENs = SDDS_GetColumnInDoubles(&SDDSin, "XmfpHENs")) ||
      !(materialData->XmfpGamD = SDDS_GetColumnInDoubles(&SDDSin, "XmfpGamD")) ||
      !(materialData->XmfpGamI = SDDS_GetColumnInDoubles(&SDDSin, "XmfpGamI"))) {
    SDDS_SetError("Some data missing from material data file");
    SDDS_PrintErrors(stdout, SDDS_VERBOSE_PrintErrors|SDDS_EXIT_PrintErrors);
  }

  printf("****************\n");
  printf("Data from file %s\n", dataFile);
  printf("%10s %4s %7s %7s %7s %7s %7s %7s %7s %7s %7s\n",
         "Material", "Z", "A", "rho", "X0", "R0", "GRNs", "MIDs", "HENs", 
         "GamD", "GamI");
  for (i=0; i<materialData->nMaterials; i++) {
    printf("%10s %4ld %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
           materialData->material[i],
           (long)materialData->Z[i],
           materialData->A[i],
           materialData->rho[i],
           materialData->X0[i],
           materialData->R0[i],
           materialData->XmfpGRNs[i],
           materialData->XmfpMIDs[i],
           materialData->XmfpHENs[i],
           materialData->XmfpGamD[i],
           materialData->XmfpGamI[i]);
  } 
  printf("****************\n");
 
  if (SDDS_ReadPage(&SDDSin)>0) 
    bomb("Material data file has multiple pages, which is not allowed.", NULL);
  SDDS_Terminate(&SDDSin);
  
  return 1;
}

