# to run PowerShell scipt, type: .\script_name on PowerShell terminal

$INPUTFILE  = "6x2mm_cartilage_0_1x_0_15mm_biased_mesh_384elem_Cp460_Creep"
$JOBNAME   	= $INPUTFILE
$SRC        = "../src/uel_pegel.for"
$NPROC      = 1

clear

echo "ABAQUS JOB RUNNING: $JOBNAME"
echo "UEL SUBROUTINE: $SRC"

abaqus interactive double analysis ask_delete=off input=$INPUTFILE job=$JOBNAME user=$SRC cpus=$NPROC

echo "ABAQUS JOB $JOBNAME IS COMPLETED"