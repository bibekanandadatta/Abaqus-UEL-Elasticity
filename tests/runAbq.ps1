# to run PowerShell scipt, type: .\script_name on PowerShell terminal

$NPROC      = 1
$JOBNAME1   = "single_elem_uniaxial_U7_CPE4_elastic"
$INPUTFILE1 = "$JOBNAME1.inp"

$UEL        = "../src/uel_mech.for"

clear

echo "ABAQUS JOB RUNNING: $JOBNAME1"
echo "UEL SUBROUTINE: $UEL"
abaqus interactive double analysis ask_delete=off job=$JOBNAME1 input=$INPUTFILE1 user=$UEL cpus=$NPROC
echo "ABAQUS JOB $JOBNAME1 IS COMPLETED"
