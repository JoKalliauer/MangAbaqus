# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2020 replay file
# Internal Version: 2019_09_13-17.49.31 163176
# Run by root on Wed Feb 24 06:29:57 2021
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=138.696258544922, 
    height=143.794448852539)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
o1 = session.openOdb(
    name='/home/jkalliau/Abaqus/MangAbaqus/AbaqusRuns/TL_arch3D-B32-20-loadfac-1-eps0.5.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Model: /home/jkalliau/Abaqus/MangAbaqus/AbaqusRuns/TL_arch3D-B32-20-loadfac-1-eps0.5.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       4
#: Number of Node Sets:          4
#: Number of Steps:              47
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=41, frame=7 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=41, frame=7 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=41, frame=6 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=41, frame=5 )
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Step-42', frame=7)
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Step-18', frame=1)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='U', outputPosition=NODAL, refinement=(INVARIANT, 
    'Magnitude'), )
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U1'), )
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Step-42', frame=7)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'), )
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Step-42', frame=7)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U3'), )
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Step-42', frame=7)
session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
    variableLabel='U', outputPosition=NODAL, refinement=(INVARIANT, 
    'Magnitude'), )
session.viewports['Viewport: 1'].odbDisplay.display.setValues(
    plotState=SYMBOLS_ON_DEF)
session.viewports['Viewport: 1'].odbDisplay.setSymbolVariable(
    variableLabel='UR', outputPosition=NODAL, vectorQuantity=RESULTANT, )
session.viewports[session.currentViewportName].odbDisplay.setFrame(
    step='Step-42', frame=7)
session.viewports['Viewport: 1'].odbDisplay.setSymbolVariable(
    variableLabel='U', outputPosition=NODAL, vectorQuantity=RESULTANT, )
session.viewports['Viewport: 1'].odbDisplay.setSymbolVariable(
    variableLabel='SM', outputPosition=INTEGRATION_POINT, 
    tensorQuantity=ALL_DIRECT_COMPONENTS, )
session.viewports['Viewport: 1'].view.setValues(nearPlane=14.5017, 
    farPlane=25.0207, width=10.3897, height=5.47435, cameraPosition=(14.5169, 
    12.1973, 11.5545), cameraUpVector=(-0.654842, 0.571747, -0.494255), 
    cameraTarget=(3.16067, 0.841038, 0.198282))
session.viewports['Viewport: 1'].view.setValues(nearPlane=15.6363, 
    farPlane=23.6889, width=11.2026, height=5.90267, cameraPosition=(2.41759, 
    15.6536, 13.1472), cameraUpVector=(-0.535406, 0.240228, -0.80971), 
    cameraTarget=(3.10457, 0.857063, 0.205667))
session.viewports['Viewport: 1'].view.setValues(nearPlane=14.9778, 
    farPlane=24.4057, width=10.7308, height=5.65408, cameraPosition=(-0.299188, 
    13.2041, 15.1339), cameraUpVector=(-0.255182, 0.43979, -0.861085), 
    cameraTarget=(3.10554, 0.857934, 0.204961))
session.viewports['Viewport: 1'].view.setValues(nearPlane=15.6782, 
    farPlane=23.5756, width=11.2326, height=5.91848, cameraPosition=(1.65031, 
    17.0628, 11.2601), cameraUpVector=(-0.0723867, 0.248226, -0.965994), 
    cameraTarget=(3.10773, 0.862274, 0.200604))
session.viewports['Viewport: 1'].view.setValues(nearPlane=14.8007, 
    farPlane=24.8403, width=10.6039, height=5.58721, cameraPosition=(7.14509, 
    6.23293, 18.6751), cameraUpVector=(-0.248821, 0.808174, -0.533801), 
    cameraTarget=(3.09577, 0.885841, 0.184468))
session.viewports['Viewport: 1'].view.setValues(nearPlane=14.4637, 
    farPlane=25.095, width=10.3625, height=5.45999, cameraPosition=(13.4116, 
    11.2803, 13.3427), cameraUpVector=(-0.266931, 0.596382, -0.757018), 
    cameraTarget=(3.14349, 0.924277, 0.143862))
session.viewports['Viewport: 1'].view.setValues(nearPlane=14.9446, 
    farPlane=24.5096, width=10.707, height=5.64151, cameraPosition=(10.6316, 
    15.0984, 11.548), cameraUpVector=(-0.33883, 0.402297, -0.850501), 
    cameraTarget=(3.12806, 0.945475, 0.133898))
