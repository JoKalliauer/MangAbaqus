from driverConstants import *
from driverStandardMPI import StandardMPIAnalysis
import driverUtils, sys
options = {
    'SIMExt':'.sim',
    'ams':OFF,
    'analysisType':STANDARD,
    'applicationName':'analysis',
    'aqua':OFF,
    'beamSectGen':OFF,
    'biorid':OFF,
    'cavityTypes':[],
    'cavparallel':OFF,
    'complexFrequency':OFF,
    'contact':OFF,
    'cosimulation':OFF,
    'coupledProcedure':OFF,
    'cpus':2,
    'cse':OFF,
    'cyclicSymmetryModel':OFF,
    'directCyclic':OFF,
    'direct_solver':DMP,
    'dsa':OFF,
    'dynStepSenseAdj':OFF,
    'dynamic':OFF,
    'filPrt':[],
    'fils':[],
    'finitesliding':OFF,
    'foundation':OFF,
    'geostatic':OFF,
    'geotech':OFF,
    'heatTransfer':OFF,
    'impJobExpVars':{},
    'importJobList':[],
    'importer':OFF,
    'importerParts':OFF,
    'includes':['eccenCompressionBeam-B32OSH-10-len-5-ecc-0.1-loadfac-1-model.inp'],
    'initialConditionsFile':OFF,
    'input':'eccenCompressionBeam-B32OSH-10-len-5-ecc-0.1-loadfac-1',
    'inputFormat':INP,
    'interactive':None,
    'job':'eccenCompressionBeam-B32OSH-10-len-5-ecc-0.1-loadfac-1',
    'keyword_licenses':[],
    'lanczos':OFF,
    'libs':[],
    'magnetostatic':OFF,
    'massDiffusion':OFF,
    'modifiedTet':OFF,
    'moldflowFiles':[],
    'moldflowMaterial':OFF,
    'mp_file_system':(DETECT, DETECT),
    'mp_head_node':('lws62.imws.tuwien.ac.at', 'lws62', '128.131.98.62', 'jkalliau-z87m-d3h', '192.168.123.1'),
    'mp_host_list':(('jkalliau-z87m-d3h', 2),),
    'mp_mode':MPI,
    'mp_mode_requested':MPI,
    'mp_mpi_validate':OFF,
    'mp_mpirun_path':'/opt/abaqus/SimulationServices/V6R2019x/linux_a64/code/bin/SMAExternal/pmpi/bin/mpirun',
    'mp_rsh_command':'ssh -n -l jkalliau %H %C',
    'multiphysics':OFF,
    'noDmpDirect':[],
    'noMultiHost':[],
    'noMultiHostElemLoop':['matrixgenerate'],
    'no_domain_check':1,
    'outputKeywords':ON,
    'parameterized':OFF,
    'partsAndAssemblies':ON,
    'parval':OFF,
    'postOutput':OFF,
    'preDecomposition':ON,
    'restart':OFF,
    'restartEndStep':OFF,
    'restartIncrement':0,
    'restartStep':0,
    'restartWrite':OFF,
    'rezone':OFF,
    'runCalculator':OFF,
    'soils':OFF,
    'soliter':OFF,
    'solverTypes':['DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT', 'DIRECT'],
    'standard_parallel':ALL,
    'staticNonlinear':ON,
    'steadyStateTransport':OFF,
    'step':ON,
    'stepSenseAdj':OFF,
    'subGen':OFF,
    'subGenLibs':[],
    'subGenTypes':[],
    'submodel':OFF,
    'substrLibDefs':OFF,
    'substructure':OFF,
    'symmetricModelGeneration':OFF,
    'thermal':OFF,
    'tmpdir':'/tmp',
    'tracer':OFF,
    'unsymm':OFF,
    'visco':OFF,
    'xplSelect':OFF,
}
analysis = StandardMPIAnalysis(options)
status = analysis.run()
sys.exit(status)
