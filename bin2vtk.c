#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>

#define LONGSTRING 200


char FilePrefix[LONGSTRING], MeshFilePrefix[LONGSTRING];
int Suffix_TMin, Suffix_TMax, Suffix_TDelta, NumFiles, DataType, ND, IntFormat = 0, ScalarTag = 0, VectorTag = 0;

void ProcessCommandLineArguements(int argc, char *argv[]);
void PrintUsage(void);
void ProcessTracerData(void);
void ProcessUnstructuredData(void);
void ReadCoordinateData(int *NumNodes, double **Coordinates);
void ReadConnectivityData(int *NumElements, int **Connectivity);
void ProcessCartesianData(void);
void FatalError(char *text, ...);

int main(int argc, char *argv[]) {
	
    ProcessCommandLineArguements(argc, argv);
	
    if(DataType == 0) 
        ProcessTracerData();
    else if(DataType >= 1 && DataType <= 4) 
        ProcessUnstructuredData();
    else if(DataType >= 5 && DataType <= 8)
        ProcessCartesianData();   
	
    return 0;
	
}

void ProcessCommandLineArguements(int argc, char *argv[]) {  
	
    if(argc < 2) 
        PrintUsage();
    else if(argc < 7 || argc > 9) 
        FatalError("Unsupported number of command line arguements\n");
    
    /* Parse argv for parameters */
    DataType = atoi(argv[1]);
    ND = atoi(argv[2]);
    sprintf(FilePrefix, "%s", argv[3]);
    Suffix_TMin = atoi(argv[4]);
    Suffix_TMax = atoi(argv[5]);
    Suffix_TDelta = atoi(argv[6]);
    if(Suffix_TMax != Suffix_TMin) {
        if(Suffix_TDelta < 1)
            FatalError("Suffix_TDelta < 1");
        NumFiles = (Suffix_TMax - Suffix_TMin) / Suffix_TDelta + 1;
    }
    else
        NumFiles = 1;
    /* Determine mesh file */	
    if(DataType != 0)
        sprintf(MeshFilePrefix, "%s", argv[7]);
    /* See if integer flag set */
    if(argc == 9) {
        if(!strcmp("-i", argv[8]))  
            IntFormat = 1;
        else
            FatalError("Unrecognized flag %s\n", argv[8]);
    }
    /* Check if tracers have scalar attached to them */
    if(DataType == 0 && argc == 8) {
        if(!strcmp("-s", argv[7])) { 
            ScalarTag = 1;
            printf("ScalarTag being set to 1\n");
        }
        else if (!strcmp("-v", argv[7])) { 
            VectorTag = 1;
            printf("VectorTag being set to 1\n");
        }
        else
            FatalError("Unrecognized flag %s (Note: tracers cannot have scalar AND vector data)\n", argv[7]);
    }   
	
    if(ND != 2 && ND != 3)
        FatalError("Invalid ND");
    if(DataType < 0 && DataType > 8)
        FatalError("Invalid DataType");
	
}

void PrintUsage(void) {
	
    printf("\nUsage:\n");
    printf("  bin2vtk DataType ND FilePrefix Start End Delta (MeshFilePrefix) (-flags)\n\n");
    printf("Description:\n");
    printf("  Converts FilePrefix.#.bin to FilePrefix.#.vtk,\n");
    printf("  where # varies from Start to End in increments of Delta.\n");   
    printf("  ND should 2 or 3 depending if data is 2D or 3D\n");
    printf("  flags:\n");
    printf("  -i: used for integer valued data\n");
    printf("  -s: used if tracer data has scalar assigned to each tracer\n");
    printf("  -v: used if tracer data has vector assigned to each tracer\n\n");
	
    printf("  Supported values for DataType:\n");
    printf("  0: Tracer position data\n");
    printf("  1: Scalar unstructured node data*\n");
    printf("  2: Vector unstructured node data*\n");
    printf("  3: Scalar unstructured element data*\n");
    printf("  4: Vector unstructured element data*\n");
    printf("  5: Scalar Cartesian node data**\n");
    printf("  6: Vector Cartesian node data**\n");
    printf("  7: Scalar Cartesian element data**\n");
    printf("  8: Vector Cartesian element data**\n");
	
    printf("  *Requires files MeshFilePrefix_coordinates.bin and MeshFilePrefix_connectivity.bin\n");
    printf("  **Requires file MeshFilePrefix_Cartesian.bin\n");
	
    exit(1);
	
}


void ProcessTracerData(void) {
	
    int i, j, num_tracers;
    double X[6], *ScalarArray, *Vx, *Vy, *Vz, time;
    char InFile[LONGSTRING], OutFile[LONGSTRING];
    FILE *InFileID, *OutFileID;
	
    for(i = 0; i < NumFiles; i++) {
        /* File names */
        sprintf(InFile, "%s.%d.bin", FilePrefix, Suffix_TMin + i * Suffix_TDelta);
        sprintf(OutFile, "%s.%d.vtk", FilePrefix, Suffix_TMin + i * Suffix_TDelta);
		
        printf("Converting %s to %s...", InFile, OutFile);
        fflush(stdout);
		
        /* Open files */
        if((InFileID = fopen(InFile, "rb")) == NULL) 
            FatalError("Could not open %s", InFile);
        if((OutFileID = fopen(OutFile, "w")) == NULL) 
            FatalError("Could not open %s", OutFile);
		
        /* Read time stamp */
        if(fread(&time, sizeof(double), 1, InFileID) < 1)
            FatalError("Could not read time stamp from file %s", InFile);      
		
        /* Read data to determine how many points there are */
        num_tracers = 0;
        if(ScalarTag)
            while(!(fread(X, sizeof(double), 4, InFileID) < 4))
                num_tracers++;
        else if(VectorTag)
            while(!(fread(X, sizeof(double), 6, InFileID) < 6))
                num_tracers++;
        else
            while(!(fread(X, sizeof(double), 3, InFileID) < 3))
                num_tracers++;
        
        printf("Number of tracers: %d\n", num_tracers);
        fflush(stdout);
        
        if(fseek(InFileID, sizeof(double), SEEK_SET))
            FatalError("Could not reset stream to start of coordinate values");
		
        /* Write VTK header */
        fprintf(OutFileID, "# vtk DataFile Version 3.0\n");
        fprintf(OutFileID, "t = %g\n", time);
        fprintf(OutFileID, "ASCII\n");
        fprintf(OutFileID, "DATASET POLYDATA\n");
        fprintf(OutFileID, "\nPOINTS %d float\n", num_tracers);
		
        /* Allocated memeory for attributes as needed */
        if(ScalarTag) {
            if((ScalarArray = (double *)calloc(num_tracers, sizeof(double))) == NULL)
                FatalError("calloc failed for ScalarArray in function ProcessTracerData()");
        }
        else if(VectorTag) {
            if((Vx = (double *)calloc(num_tracers, sizeof(double))) == NULL)
                FatalError("calloc failed for Vx in function ProcessTracerData()");
            if((Vy = (double *)calloc(num_tracers, sizeof(double))) == NULL)
                FatalError("calloc failed for Vy in function ProcessTracerData()");
            if((Vz = (double *)calloc(num_tracers, sizeof(double))) == NULL)
                FatalError("calloc failed for Vz in function ProcessTracerData()");
		}
        /* Read/Write data */
        for(j = 0; j < num_tracers; j++) {
            if(fread(X, sizeof(double), 3, InFileID) < 3) 
                FatalError("Could not read location of tracer %d of %d from %s", j + 1, num_tracers, InFile);
            if(ScalarTag) {
                if(fread(&ScalarArray[j], sizeof(double), 1, InFileID) < 1)
                    FatalError("Could not read scalar value for tracer %d of %d from %s", j + 1, num_tracers, InFile);
            }
            else if(VectorTag) {
                if(fread(&Vx[j], sizeof(double), 1, InFileID) < 1)
                    FatalError("Could not read vector value 1 for tracer %d of %d from %s", j + 1, num_tracers, InFile);
                if(fread(&Vy[j], sizeof(double), 1, InFileID) < 1)
                    FatalError("Could not read vector value 2 for tracer %d of %d from %s", j + 1, num_tracers, InFile);
                if(fread(&Vz[j], sizeof(double), 1, InFileID) < 1)
                    FatalError("Could not read vector value 3 for tracer %d of %d from %s", j + 1, num_tracers, InFile);
            }
            if(ND == 3)
                fprintf(OutFileID, "%.9f %.9f %.9f\n", X[0], X[1], X[2]);
            else
                fprintf(OutFileID, "%.9f %.9f 0.0\n", X[0], X[1]);						  
        }
		
        if(ScalarTag) {
            fprintf(OutFileID, "\nPOINT_DATA %d\n", num_tracers);
            fprintf(OutFileID, "SCALARS scalar float 1\n");
            fprintf(OutFileID, "LOOKUP_TABLE default\n");
            for(j = 0; j < num_tracers; j++) 
                fprintf(OutFileID, "%.9f ", ScalarArray[j]);
        }
        else if(VectorTag) {
            fprintf(OutFileID, "\nPOINT_DATA %d\n", num_tracers);
            fprintf(OutFileID, "VECTORS vector float\n");
            /*fprintf(OutFileID, "LOOKUP_TABLE default\n");*/
            for(j = 0; j < num_tracers; j++) 
                fprintf(OutFileID, "%.9f %.9f %.9f\n", Vx[j], Vy[j], Vz[j]);
        }
		
		/* Free memeory for attributes as needed */
        if(ScalarTag) 
            free(ScalarArray);  
        else if(VectorTag) {
            free(Vx);
            free(Vy);
            free(Vz);
        }
        
        printf("OK!\n");
        fflush(stdout);
		
        /* Close files */
        fclose(InFileID);  
        fclose(OutFileID);
    }
	
}

void ProcessUnstructuredData(void) {
	
    int i, j, NumNodes, NumElements, **Connectivity, ISV, IVV[3];
    double **Coordinates, FSV, FVV[3], time;
    char InFile[LONGSTRING], OutFile[LONGSTRING];
    FILE *InFileID, *OutFileID;
	
    /* Open coordinate data file */
    sprintf(InFile, "%s_coordinates.bin", MeshFilePrefix);
    if((InFileID = fopen(InFile, "rb")) == NULL) 
        FatalError("Could not open %s", InFile);
	
    /* Read number of nodes */
    if(fread(&NumNodes, sizeof(int), 1, InFileID) < 1)
        FatalError("Could not read number of nodes from %s", InFile);
	
    printf("Reading coordinates for %d nodes...", NumNodes);
    fflush(stdout);
	
    /* Allocate memory for Coordinates array */
    if((Coordinates = (double **)malloc(NumNodes * sizeof(double *))) == NULL)
        FatalError("Malloc failed for Coordinates");
    for(i = 0; i < NumNodes; i++)
        if((Coordinates[i] = (double *)malloc(3 * sizeof(double))) == NULL)
            FatalError("Malloc failed for Coordinates[%d]", i);
	
    /* Read coordinates */
    for(i = 0; i < NumNodes; i++) /* { */ 
        if(fread(Coordinates[i], sizeof(double), 3, InFileID) < 3)
            FatalError("Could read Coordinates[%d] from %s", i, InFile);
	
    printf("OK!\n");
	
    /* Close coordinates file */
    fclose(InFileID);
	
    /* Open connectivity data file */
    sprintf(InFile, "%s_connectivity.bin", MeshFilePrefix);
    if((InFileID = fopen(InFile, "rb")) == NULL) 
        FatalError("Could not open %s", InFile);
	
    /* Read number of elements */
    if(fread(&NumElements, sizeof(int), 1, InFileID) < 1)
        FatalError("Could not read number of elements from %s", InFile);
	
    printf("Reading connectivity for %d elements...", NumElements);
    fflush(stdout);
	
    /* Allocate memory for Connectivity array */
    if((Connectivity = (int **)malloc(NumElements * sizeof(int *))) == NULL)
        FatalError("Malloc failed for Connectivity");
    for(i = 0; i < NumElements; i++)
        if((Connectivity[i] = (int *)malloc(4 * sizeof(int))) == NULL)
            FatalError("Malloc failed for Connectivity[%d]", i);
	
    /* Read connectivity */
    for(i = 0; i < NumElements; i++) /* { */
        if(fread(Connectivity[i], sizeof(int), 4, InFileID) < 4)
            FatalError("Could not read Connectivity[%d] from %s", i, InFile);
	
    printf("OK!\n");
	
    /* Close connectivity file */
    fclose(InFileID);
	
    /* Process field data files */
    for(i = 0; i < NumFiles; i++) {
        /* File names */
        sprintf(InFile, "%s.%d.bin", FilePrefix, Suffix_TMin + i * Suffix_TDelta);
        sprintf(OutFile, "%s.%d.vtk", FilePrefix, Suffix_TMin + i * Suffix_TDelta);
		
        printf("Converting %s to %s...\n", InFile, OutFile);
        fflush(stdout);
		
        /* Open files */
        if((InFileID = fopen(InFile, "rb")) == NULL) 
            FatalError("Could not open %s", InFile);
        if((OutFileID = fopen(OutFile, "w")) == NULL) 
            FatalError("Could not open %s", OutFile);
		
        /* Read time stamp */
        if(fread(&time, sizeof(double), 1, InFileID) < 1)
            FatalError("Could not read time stamp from file %s", InFile); 
		
        printf("  t = %f\n", time);
        fflush(stdout);
		
        /* Write VTK header */
        fprintf(OutFileID, "# vtk DataFile Version 3.0\n");
        fprintf(OutFileID, "t = %.9f\n", time);
        fprintf(OutFileID, "ASCII\n");
        fprintf(OutFileID, "DATASET UNSTRUCTURED_GRID\n\n");
		
        printf("  Printing coordinates...\n");
        fflush(stdout);
		
        /* Print node coordinates */
        fprintf(OutFileID, "POINTS %d float\n", NumNodes);
        for(j = 0; j < NumNodes; j++) {
            if(ND == 3)
                fprintf(OutFileID, "%.9f %.9f %.9f\n", Coordinates[j][0],  Coordinates[j][1], Coordinates[j][2]); 
            else
                fprintf(OutFileID, "%.9f %.9f 0.0\n", Coordinates[j][0],  Coordinates[j][1]); 
        }
		
        printf("  Printing connectivity...\n");
        fflush(stdout);
		
        /* Print connectivity */
        if(ND  == 3) {
            fprintf(OutFileID, "\nCELLS %d %d\n", NumElements, 5*NumElements);
            for(j = 0; j < NumElements; j++) 
                fprintf(OutFileID, "4 %d %d %d %d\n", Connectivity[j][0], Connectivity[j][1], Connectivity[j][2], Connectivity[j][3]);
        }
        else {
            fprintf(OutFileID, "\nCELLS %d %d\n", NumElements, 4*NumElements);
            for(j = 0; j < NumElements; j++) 
                fprintf(OutFileID, "3 %d %d %d\n", Connectivity[j][0], Connectivity[j][1], Connectivity[j][2]);
        }
		
        /* Print cell type (10 = VTK_TETRA, 5 = VTK_TRIANGLE) */
        fprintf(OutFileID, "\nCELL_TYPES %d\n", NumElements);
        if(ND  == 3) 
            for(j = 0; j < NumElements; j++) 
                fprintf(OutFileID, "10\n"); 
        else
            for(j = 0; j < NumElements; j++) 
                fprintf(OutFileID, "5\n");
		
        /* Print header for field data */
        if(DataType == 1 || DataType == 2)  /* Node data */
            fprintf(OutFileID, "\nPOINT_DATA %d\n", NumNodes);
        else /* Element data */
            fprintf(OutFileID, "\nCELL_DATA %d\n", NumElements);
        if(DataType == 1 || DataType == 3) { /* Scalar values */
            if(IntFormat) 	
                fprintf(OutFileID, "SCALARS scalar_field int\n");
            else /* Floating point data */
                fprintf(OutFileID, "SCALARS scalar_field float\n");
            fprintf(OutFileID, "LOOKUP_TABLE default\n");
        }
        else { /* Vector values */
            if(IntFormat) 	
                fprintf(OutFileID, "VECTORS vector_field int\n");
            else /* Floating point data */
                fprintf(OutFileID, "VECTORS vector_field float\n");
        }
		
		
        printf("  Printing field data...\n");
        fflush(stdout);
		
        /* Read/Write data */
        for(j = 0; j < ((DataType == 1 || DataType == 2)? NumNodes:NumElements); j++) { 
            if(DataType == 1 || DataType == 3) { /* Scalar data */
                if(IntFormat) {
                    if(fread(&ISV, sizeof(int), 1, InFileID) < 1)
                        FatalError("Could not read integer scalar value %d of %d from file %s", j + 1, (DataType == 1 || DataType == 2)? NumNodes:NumElements, InFile);
                    fprintf(OutFileID, "%d\n", ISV);
                }
                else { /* Floating point data */
                    if(fread(&FSV, sizeof(double), 1, InFileID) < 1)
                        FatalError("Could not read float scalar value %d of %d from file %s", j + 1, (DataType == 1 || DataType == 2)? NumNodes:NumElements, InFile);
                    fprintf(OutFileID, "%.9f\n", FSV);
                }
            }
            else { /* Vector data */
                if(IntFormat) {
                    if(fread(IVV, sizeof(int), 3, InFileID) < 3)
                        FatalError("Could not read integer vector value %d of %d from file %s", j + 1, (DataType == 1 || DataType == 2)? NumNodes:NumElements, InFile);
                    if(ND == 3)
                        fprintf(OutFileID, "%d %d %d\n", IVV[0], IVV[1], IVV[2]);
                    else /* ND == 2 */
                        fprintf(OutFileID, "%d %d\n", IVV[0], IVV[1]);
                }
                else { /* Floating point data */
                    if(fread(&FVV, sizeof(double), 3, InFileID) < 3)
                        FatalError("Could not read float vector value %d of %d from file %s", j + 1, (DataType == 1 || DataType == 2)? NumNodes:NumElements, InFile);
                    if(ND == 3)
                        fprintf(OutFileID, "%.9f %.9f %.9f\n", FVV[0], FVV[1], FVV[2]);
                    else /* ND == 2 */
                        fprintf(OutFileID, "%.9f %.9f 0.0\n", FVV[0], FVV[1]);
                }
            }
        }
		
        printf("OK!\n");
        fflush(stdout);
		
        /* Close files */
        fclose(InFileID);
        fclose(OutFileID); 
    } 
	
    /* Free memory */
    for(i = 0; i < NumNodes; i++) 
        free(Coordinates[i]);
    free(Coordinates);
    for(i = 0; i < NumElements; i++)   
        free(Connectivity[i]);
    free(Connectivity);
	
}


void ProcessCartesianData(void) {
	
    int ii, i, j, k, imax, jmax, kmax, XRes, YRes, ZRes, NumNodes, NumElements;
    double time, SV, VV[3], XMin, XMax, XDelta, YMin, YMax, YDelta, ZMin, ZMax, ZDelta;
    char InFile[LONGSTRING], OutFile[LONGSTRING];
    FILE *InFileID, *OutFileID;
	
    /* Open Cartesian mesh information file */
    sprintf(InFile, "%s_Cartesian.bin", MeshFilePrefix);
    if((InFileID = fopen(InFile, "rb")) == NULL) 
        FatalError("Could not open file %s", InFile);
	
    /* Read mesh parameters */
    if(fread(&XMin, sizeof(double), 1, InFileID) < 1 ||
       fread(&XMax, sizeof(double), 1, InFileID) < 1 ||
       fread(&XRes, sizeof(int),    1, InFileID) < 1 ||
       fread(&YMin, sizeof(double), 1, InFileID) < 1 ||
       fread(&YMax, sizeof(double), 1, InFileID) < 1 ||
       fread(&YRes, sizeof(int),    1, InFileID) < 1 ||
       fread(&ZMin, sizeof(double), 1, InFileID) < 1 ||
       fread(&ZMax, sizeof(double), 1, InFileID) < 1 ||
       fread(&ZRes, sizeof(int),    1, InFileID) < 1)
        FatalError("Could not read Catesian mesh parameters from %s", InFile);
    /* Close file */
    fclose(InFileID);
    
    /* Derived parameters */
    XDelta = (XMax - XMin) / (XRes - 1);
    YDelta = (YMax - YMin) / (YRes - 1);
    if(ND == 3)
        ZDelta = (ZMax - ZMin) / (ZRes - 1);
    else {
        ZDelta = 0;
        ZRes = 1;
    }
    NumNodes = XRes * YRes * ZRes;
    if(ND == 3)
        NumElements = (XRes - 1) * (YRes - 1) * (ZRes - 1);
    else
        NumElements = (XRes - 1) * (YRes - 1);
	
    for(ii = 0; ii < NumFiles; ii++) {
        /* File names */
        sprintf(InFile, "%s.%d.bin", FilePrefix, Suffix_TMin + ii * Suffix_TDelta);
        sprintf(OutFile, "%s.%d.vtk", FilePrefix, Suffix_TMin + ii * Suffix_TDelta);
        
        printf("Converting %s to %s...", InFile, OutFile);
        fflush(stdout);
		
        /* Open files */
        if((InFileID = fopen(InFile, "rb")) == NULL) 
            FatalError("Could not open %s", InFile);
        if((OutFileID = fopen(OutFile, "w")) == NULL) 
            FatalError("Could not open %s", OutFile);
		
        /* Read time stamp */
        if(fread(&time, sizeof(double), 1, InFileID) < 1)
            FatalError("Could not read time stamp from file %s", InFile);     
		
        /* Write VTK header */
        fprintf(OutFileID, "# vtk DataFile Version 3.0\n");
        fprintf(OutFileID, "t = %.9f\n", time);
        fprintf(OutFileID, "ASCII\n");
        fprintf(OutFileID, "DATASET STRUCTURED_POINTS\n");
        fprintf(OutFileID, "DIMENSIONS %d %d %d\n", XRes, YRes, ZRes);
        fprintf(OutFileID, "ORIGIN %f %f %f\n", XMin, YMin, ZMin);
        fprintf(OutFileID, "SPACING %.9f %.9f %.9f\n", XDelta, YDelta, ZDelta);
        
        /* Print header for field data */
        if(DataType == 5 || DataType == 6)  /* Node data */
            fprintf(OutFileID, "\nPOINT_DATA %d\n", NumNodes);
        else if(DataType == 7 || DataType == 8) /* Element data */
            fprintf(OutFileID, "\nCELL_DATA %d\n", NumElements);
        else
            FatalError("Unexpected value, DataType = %d");
        if(DataType == 5 || DataType == 7) { /* Scalar values */
            if(IntFormat) 	
                fprintf(OutFileID, "SCALARS scalar_field int\n");
            else /* Floating point data */
                fprintf(OutFileID, "SCALARS scalar_field float\n");
            fprintf(OutFileID, "LOOKUP_TABLE default\n");
        }
        else { /* Vector values */
            if(IntFormat) 	
                fprintf(OutFileID, "VECTORS vector_field int\n");
            else /* Floating point data */
                fprintf(OutFileID, "VECTORS vector_field float\n");
        }
        
        /* Read/Write data */
        imax = (DataType == 5 || DataType == 6)?XRes:(XRes-1);
        jmax = (DataType == 5 || DataType == 6)?YRes:(YRes-1);
        if(ND == 3)
            kmax = (DataType == 5 || DataType == 6)?ZRes:(ZRes-1);
        else
            kmax = 1;
        for(k = 0; k < kmax; k++) 
            for(j = 0; j < jmax; j++) 
                for(i = 0; i < imax; i++) 
                    if(DataType == 5 || DataType == 7) { /* Scalar data */
                        if(fread(&SV, sizeof(double), 1, InFileID) < 1)
                            FatalError("Could read scalar value (%d, %d, %d) from file %s, SV = %f", i, j, k, InFile, SV);
                        fprintf(OutFileID, "%.9f\n", SV);
                    }
                    else { /* Vector data */
                        if(fread(VV, sizeof(double), 3, InFileID) < 3)
                            FatalError("Could not read vector value (%d, %d, %d) from file %s", i, j, k, InFile);
                        fprintf(OutFileID, "%.9f %.9f %.9f\n", VV[0], VV[1], VV[2]); 
                    }
        
        printf("OK!\n");
        fflush(stdout);
        
        /* Close files */
        fclose(InFileID);
        fclose(OutFileID);    
    }
	
}


void FatalError(char *text, ...) {
	va_list ap;
	char *p, *sval;
	int ival;
	double dval;
	
	fprintf(stderr, "\nERROR:\n");
	va_start(ap, text);
	for(p = text; *p; p++) {
		if(*p != '%') {
			putc(*p, stderr);
			continue;
		}
		switch (*++p) {
			case 'd':
				ival = va_arg(ap, int);
				fprintf(stderr, "%d", ival);
				break;
			case 'f':
				dval = va_arg(ap, double);
				fprintf(stderr, "%f", dval);
				break;
			case 's':
				for(sval = va_arg(ap, char *); *sval; sval++)
					putc(*sval, stderr);
				break;
			default:
				putc(*p, stderr);
				break;
		}
	}
	va_end(ap);
	
	fprintf(stderr, "\n");
	fflush(stderr);
	
	exit(1);
}
