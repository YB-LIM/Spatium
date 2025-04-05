from abaqusGui import *
from abaqusConstants import ALL
import osutils, os


###########################################################################
# Class definition
###########################################################################

class Spatium_v2_plugin(AFXForm):

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, owner):
        
        # Construct the base class.
        #
        AFXForm.__init__(self, owner)
        self.radioButtonGroups = {}

        self.cmd = AFXGuiCommand(mode=self, method='excecute',
            objectName='Spatium_Kernel', registerQuery=False)
        pickedDefault = ''
        self.Is_PorousKw = AFXStringKeyword(self.cmd, 'Is_Porous', True, 'Composite')
        self.LKw = AFXFloatKeyword(self.cmd, 'L', True, 0.2)
        self.L_meshKw = AFXFloatKeyword(self.cmd, 'L_mesh', True, 0.008)
        self.r_avgKw = AFXFloatKeyword(self.cmd, 'r_avg', True, 0.04)
        self.r_stdKw = AFXFloatKeyword(self.cmd, 'r_std', True, 0.005)
        self.VoF_tarKw = AFXFloatKeyword(self.cmd, 'VoF_tar', True, 0.2)
        self.min_distanceKw = AFXFloatKeyword(self.cmd, 'min_distance', True, 0.001)
        self.max_iterationsKw = AFXFloatKeyword(self.cmd, 'max_iterations', True, 100000)
        self.ModeKw = AFXStringKeyword(self.cmd, 'Mode', True, 'Uniaxial Tension')
        self.DispKw = AFXFloatKeyword(self.cmd, 'Disp', True, 0.007)
        self.Odb_PathKw = AFXStringKeyword(self.cmd, 'Odb_Path', True, '')
        self.Output_PathKw = AFXStringKeyword(self.cmd, 'Output_Path', True, '')
        self.EngStrainKw = AFXFloatKeyword(self.cmd, 'EngStrain', True, 0.0)
        self.S_CompKw = AFXStringKeyword(self.cmd, 'S_Comp', True, 'S11')
        self.Plot_FlagKw = AFXBoolKeyword(self.cmd, 'Plot_Flag', AFXBoolKeyword.TRUE_FALSE, True, True)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def getFirstDialog(self):
        self.cmd.setKeywordValuesToDefaults()
        import spatium_v2DB
        return spatium_v2DB.Spatium_v2DB(self)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def doCustomChecks(self):

        # Try to set the appropriate radio button on. If the user did
        # not specify any buttons to be on, do nothing.
        #
        for kw1,kw2,d in self.radioButtonGroups.values():
            try:
                value = d[ kw1.getValue() ]
                kw2.setValue(value)
            except:
                pass
        return True

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def okToCancel(self):

        # No need to close the dialog when a file operation (such
        # as New or Open) or model change is executed.
        #
        return False

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Register the plug-in
#
thisPath = os.path.abspath(__file__)
thisDir = os.path.dirname(thisPath)

toolset = getAFXApp().getAFXMainWindow().getPluginToolset()
toolset.registerGuiMenuButton(
    buttonText='Spatium', 
    object=Spatium_v2_plugin(toolset),
    messageId=AFXMode.ID_ACTIVATE,
    icon=None,
    kernelInitString='import Spatium_Kernel',
    applicableModules=ALL,
    version='N/A',
    author='N/A',
    description='N/A',
    helpUrl='N/A'
)
